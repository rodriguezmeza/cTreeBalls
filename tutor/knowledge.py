"""Rebuildable SQLite document index and local, word-based scientific retrieval."""
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import unicodedata

SCHEMA = 2
STOP = set("a an and are as at be been by can could do does for from have how i if in is it me my of on or our please than that the these this to us was we what when where which why will with would you your about tell explain question value values make".split())


def terms(text):
    text = text.lower().replace("κ", " kappa ").replace("γ", " shear ")
    text = text.replace("lyα", " lya ").replace("lyman-alpha", " lya ").replace("lyman alpha", " lya ")
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    for old, new in (("two point", "2pcf"), ("two-point", "2pcf"),
                     ("three point", "3pcf"), ("three-point", "3pcf"),
                     ("kd-tree", "kdtree"), ("ball-tree", "balltree"),
                     ("edge correction", "edge-corrections")):
        text = text.replace(old, new)
    text = re.sub(r"\b(improve|improves|improvement|improvements|accelerate)\b", "speedup", text)
    words = re.findall(r"[a-z0-9]+(?:[_.-][a-z0-9]+)*", text)
    # Keep the full parameter name and its pieces: natural language and exact
    # option-name questions should both retrieve the relevant explanation.
    expanded = []
    for word in words:
        if word not in STOP:
            expanded.append(word)
        if "_" in word or "-" in word:
            expanded.extend(p for p in re.split(r"[_-]", word) if p not in STOP)
    return expanded


@dataclass
class Passage:
    source: str
    line: int
    title: str
    text: str
    kind: str = "document"
    topic: str = ""
    keywords: str = ""
    references: tuple = ()
    stale: bool = False

    @property
    def citation(self):
        return f"{self.source}:{self.line}"


def source_files(root):
    """Select project documentation, excluding outputs, vendors and backups."""
    candidates = [root/"README.md", root/"ENGINE_CAPABILITIES.md"]
    for directory in (root/"docs", root/"tutor/notes"):
        if directory.exists():
            for p in directory.rglob("*"):
                relative = p.relative_to(directory)
                if any(part.startswith(("_", ".")) or part.lower() in
                       {"man", "backup", "backups", "archive", "archives", "delete"}
                       for part in relative.parts):
                    continue
                if p.suffix.lower() in {".md", ".rst"}:
                    candidates.append(p)
    for sub in ("tests/python", "tests/slurm_scripts", "scripts", "support"):
        candidates.extend((root/sub).glob("README*"))
    profile = root/"capabilities/build_profile.json"
    if profile.exists():
        for name in json.loads(profile.read_text())["enabled_addon_directories"]:
            candidates.extend((root/"addons"/name).glob("README*"))
    candidates += [root/"capabilities/engines.json", root/"tutor/knowledge/questions.json"]
    result = []
    for p in sorted(set(candidates)):
        # A documentation symlink must not silently pull a private file into
        # the index or the optional model context.
        if p.is_file() and p.resolve().is_relative_to(root) and p.stat().st_size <= 2_000_000:
            result.append(p)
    return result


def document_passages(source, text):
    lines = text.splitlines()
    start, title, block, fenced = 1, source, [], False
    for number, line in enumerate(lines, 1):
        heading = re.match(r"^#{1,4}\s+(.+)", line) if not fenced else None
        rst_heading = (not fenced and number < len(lines) and len(line.strip()) > 1
                       and bool(re.fullmatch(r"[=~^-]{3,}", lines[number].strip())))
        if (heading or rst_heading) and block:
            yield Passage(source, start, title, "\n".join(block).strip())
            block, start = [], number
        if heading or rst_heading:
            title = heading.group(1) if heading else line.strip()
        if line.lstrip().startswith("```"):
            fenced = not fenced
        block.append(line)
        # Avoid splitting a small code example or a table in mid-line. Long
        # sections remain searchable and retain exact source line numbers.
        if not fenced and not line.strip() and sum(map(len, block)) > 2100:
            yield Passage(source, start, title, "\n".join(block).strip())
            block, start = [], number+1
    if any(s.strip() for s in block):
        yield Passage(source, start, title, "\n".join(block).strip())


class Knowledge:
    def __init__(self, root, database=None, rebuild=False):
        self.root = Path(root).resolve()
        database = database or self.root/"build/ia-tutor/knowledge.sqlite3"
        self.database = str(database)
        self.cache_notice = ""
        files = source_files(self.root)
        self.file_hashes = {str(p.relative_to(self.root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
        manifest = json.dumps(dict(schema=SCHEMA, root=str(self.root), files=self.file_hashes), sort_keys=True)
        try:
            if self.database != ":memory:":
                Path(self.database).parent.mkdir(parents=True, exist_ok=True)
            self._open_cache()
            self._read_cache(files, manifest, rebuild)
        except (OSError, sqlite3.Error):
            if hasattr(self, "connection"):
                self.connection.close()
            self.database = ":memory:"
            self.cache_notice = "Using an in-memory index; the checkout cache is not writable or usable."
            self._open_cache()
            self._read_cache(files, manifest, True)
        self.counts = [Counter(terms(p.title+" "+p.text+" "+p.keywords)) for p in self.passages]
        self.lengths = [sum(c.values()) for c in self.counts]
        self.average_length = sum(self.lengths)/max(len(self.lengths), 1)
        frequency = Counter(t for c in self.counts for t in c)
        self.idf = {t: math.log(1+(len(self.passages)-n+.5)/(n+.5)) for t, n in frequency.items()}

    def _open_cache(self):
        self.connection = sqlite3.connect(self.database, timeout=15)
        self.connection.execute("CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        self.connection.execute("CREATE TABLE IF NOT EXISTS passages (source TEXT, line INTEGER, title TEXT, text TEXT, kind TEXT, topic TEXT, keywords TEXT, refs TEXT, stale INTEGER)")

    def _read_cache(self, files, manifest, rebuild):
        cached = self.connection.execute("SELECT value FROM meta WHERE key='manifest'").fetchone()
        self.updated = rebuild or cached != (manifest,)
        if self.updated:
            passages = []
            for path in files:
                source, text = str(path.relative_to(self.root)), path.read_text(encoding="utf-8", errors="replace")
                if source == "tutor/knowledge/questions.json":
                    for item in json.loads(text)["questions"]:
                        marker = f'"id": "{item["id"]}"'
                        number = text[:text.index(marker)].count("\n")+1
                        changed = any(self.file_hashes.get(p) != digest for p, digest in item.get("source_sha256", {}).items())
                        refs = tuple(item.get("sources", []))
                        passages.append(Passage(source, number, item["title"], item["answer"], "answer", item["id"],
                            "\n".join(item["questions"]+[" ".join(item.get("keywords", []))]), refs, changed))
                elif source == "capabilities/engines.json":
                    for engine in json.loads(text)["engines"]:
                        marker = f'"name": "{engine["name"]}"'
                        number = text[:text.index(marker)].count("\n")+1
                        passages.append(Passage(source, number, engine["name"],
                            engine["geometry"]+"\n"+engine["correlations"]+"\n"+engine["usage"], "engine"))
                else:
                    passages.extend(document_passages(source, text))
            with self.connection:
                self.connection.execute("DELETE FROM passages")
                self.connection.executemany("INSERT INTO passages VALUES (?,?,?,?,?,?,?,?,?)", [
                    (p.source, p.line, p.title, p.text, p.kind, p.topic, p.keywords, json.dumps(p.references), int(p.stale)) for p in passages])
                self.connection.execute("INSERT OR REPLACE INTO meta VALUES ('manifest',?)", (manifest,))
        self.passages = [Passage(*row[:7], tuple(json.loads(row[7])), bool(row[8]))
                         for row in self.connection.execute("SELECT * FROM passages ORDER BY source,line")]

    def search(self, question, limit=5, kind=None):
        query = set(terms(question))
        ranked = []
        for passage, counts, length in zip(self.passages, self.counts, self.lengths):
            if kind and passage.kind != kind:
                continue
            overlap = query & counts.keys()
            if not overlap:
                continue
            norm = 1.2*(.25+.75*length/max(self.average_length, 1))
            score = sum(self.idf[t]*counts[t]*2.2/(counts[t]+norm) for t in overlap)
            headings = set(terms(passage.title))
            score += sum(2*self.idf[t] for t in query & headings)
            if passage.kind == "answer":
                aliases = set(terms(passage.keywords))
                score += sum(1.3*self.idf[t] for t in query & aliases)
                if any(query == set(terms(alias)) for alias in passage.keywords.splitlines()):
                    score += 50
                if passage.stale:
                    score *= .65
            coverage = sum(self.idf.get(t, 3) for t in overlap)/max(sum(self.idf.get(t, 3) for t in query), 1)
            ranked.append((score*coverage, coverage, passage))
        return sorted(ranked, key=lambda row: (-row[0], row[2].source, row[2].line))[:limit]

    def references(self, passage):
        return [passage.citation]+list(passage.references)

    def close(self):
        self.connection.close()
