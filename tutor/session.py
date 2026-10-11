"""Scientific question answering, session context, and engine guidance."""
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

from .knowledge import terms


@dataclass
class Answer:
    question: str
    text: str
    sources: list = field(default_factory=list)
    mode: str = "offline"
    topic: str = ""
    context: dict = field(default_factory=dict)
    notice: str = ""

    def as_dict(self):
        return asdict(self)


class Tutor:
    def __init__(self, knowledge, executable=None, model=None):
        self.knowledge = knowledge
        self.root = knowledge.root
        self.executable = Path(executable).resolve() if executable else self.root/"cballs"
        self.model = model
        self.context = {}
        self.history = []
        self.last = None
        self._available = None
        self.build_notice = ""
        self.catalog = json.loads((self.root/"capabilities/engines.json").read_text())["engines"]

    def set_context(self, text):
        allowed = {"field": {"kappa", "shear", "lya", "counts"},
                   "geometry": {"sphere", "flat", "3d", "1d"},
                   "statistic": {"2pcf", "3pcf", "both"},
                   "catalog": None, "machine": None, "points": None, "threads": None}
        pending = dict(self.context)
        for token in shlex.split(text):
            key, separator, value = token.partition("=")
            if not separator or key not in allowed or not value:
                raise ValueError("Use field=, geometry=, statistic=, catalog=, machine=, points= or threads=.")
            if allowed[key] is not None and value not in allowed[key]:
                raise ValueError(f"{key} must be one of: {', '.join(sorted(allowed[key]))}")
            if key in {"points", "threads"} and (not value.isdigit() or int(value) < 1):
                raise ValueError(f"{key} must be a positive integer.")
            if len(value) > 1000:
                raise ValueError("Context value is too long.")
            pending[key] = int(value) if key in {"points", "threads"} else value
        self.context = pending
        return self.context

    def infer_context(self, question):
        words = set(terms(question))
        fields = words & {"kappa", "shear", "lya", "counts"}
        if len(fields) == 1:
            self.context["field"] = fields.pop()
        for word, value in (("spherical", "sphere"), ("sphere", "sphere"),
                            ("flat", "flat"), ("cartesian", "3d"), ("1d", "1d")):
            if word in words:
                self.context["geometry"] = value
        stats = words & {"2pcf", "3pcf"}
        if stats:
            self.context["statistic"] = "both" if len(stats) == 2 else next(iter(stats))
        if "perlmutter" in words:
            self.context["machine"] = "Perlmutter"
        if words & {"des", "desy3", "des-y3", "takahashi"}:
            self.context["catalog"] = "DES/Takahashi"
            self.context.setdefault("geometry", "sphere")

    def probe(self, option):
        if option not in {"print-search-methods", "make-info"}:
            raise ValueError("Unsupported information query.")
        if not self.executable.is_file():
            raise ValueError("No cballs executable found; source capabilities are not a compiled-build check.")
        # These are fixed information commands, not commands generated from
        # questions or documents. A temporary cwd isolates legacy output paths.
        with tempfile.TemporaryDirectory(prefix="cballs-tutor-info-") as directory:
            process = subprocess.run([str(self.executable), "options="+option],
                cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10)
        if process.returncode not in (0, 1) or not process.stdout:
            raise ValueError("The cballs information query failed: "+process.stderr[-300:])
        return process.stdout

    def available(self):
        if self._available is None:
            try:
                output = self.probe("print-search-methods")
                names = set(re.findall(r"^- ([\w-]+) \(id=\d+\)", output, re.MULTILINE))
                if not names:
                    raise ValueError("The executable returned no recognizable search-method list.")
                self._available = names
                self.build_notice = "Availability checked against "+str(self.executable)
            except (OSError, ValueError, subprocess.SubprocessError) as exc:
                self._available = set()
                self.build_notice = str(exc)
        return self._available

    def engines(self, family=""):
        available = self.available()
        rows = [e for e in self.catalog if not family or family.lower() in
                (e["name"]+" "+e["family"]+" "+e["geometry"]).lower()]
        return self.build_notice+"\n"+"\n".join(
            f"  {e['name']} [{'compiled' if e['name'] in available else 'not verified in executable'}]\n    {e['correlations']}"
            for e in rows)

    def recommend(self, question="Engine selection"):
        family = self.context.get("field")
        geometry = self.context.get("geometry")
        statistic = self.context.get("statistic", "3pcf" if family in {"shear", "kappa"} else "both")
        missing = [name for name, value in (("field", family), ("geometry", geometry)) if not value]
        if missing:
            text = ("Engine choice depends on the field and geometry. Please specify "+", ".join(missing)+
                    ". For example:\n/context field=shear geometry=sphere statistic=3pcf catalog=DESY3 threads=16\n"
                    "Then ask which engine to compare. For Lyα, distinguish 3D forest pixels from radial-only 1D statistics.")
            return Answer(question, text, ["docs/search_methods.rst"], topic="engine-selection", context=dict(self.context))
        if family == "shear" and geometry != "sphere":
            return Answer(question, "The current native shear profile is spherical. Do not pass flat projected coordinates to a sphere engine. For a flat catalogue, use the private dual node comparison workflow with its flat geometry or explicitly validate a spherical coordinate conversion.",
                          ["docs/shear.rst", "docs/ACTIVE_PROFILE.md"], context=dict(self.context))
        if family == "lya" and geometry not in {"1d", "3d"}:
            return Answer(question, "Lyα forest engines need observer-centered 3D pixel positions and integer forest IDs, or a deliberately selected radial 1D estimator. A sky-only angular catalogue is insufficient.", ["docs/lyman_alpha.rst"], context=dict(self.context))
        names, detail = [], ""
        if family == "shear":
            names = [f"{tree}-shear-sphere-2balls-omp" for tree in ("balltree", "kdtree", "octree")]
            detail = ("For DES/Takahashi spherical shear, start by comparing ball-tree and KD-tree; include the octree at matching measured error. "
                      "Use the local-east-north convention when that is the catalogue's convention. For 3PCF, opt-in shear-pivot-reuse uses phase 3 and bin allowance 0.5 in the shipped profile; nsmooth defaults to 8. "
                      "A historical tighter trial with bin allowance 0.05 and nsmooth=1 on one million DES objects gave a measured 1.61x ball-tree speedup with 0.764% raw and 0.283% corrected L2 error; it is not an error guarantee for your catalogue. "
                      "Reuse does not accelerate the separate 2PCF path.")
        elif family == "kappa" and geometry in {"sphere", "flat"}:
            names = [f"{tree}-2balls-omp" for tree in ("balltree", "kdtree", "octree")]
            detail = ("Compare these scalar angular engines on the same retained catalogue, bin edges and multipoles. "
                      "At theta=1 the DES 100,000-object example had 20–28% raw 3PCF error with ordinary traversal. "
                      "Scalar reuse (phase 0.1, bin allowance 0) greatly reduced that error but was slower on that sample. "
                      "Neither reuse nor a particular tree is always faster. Check edge-corrected errors separately. "
                      "For the sphere use unit directions and chord bins; do not recenter the sky.")
        elif family == "lya":
            stem = "lya-1d-tree-" if geometry == "1d" else "lya-los-tree-"
            suffix = "2pcf-3pcf" if statistic == "both" else statistic
            names = ([stem+"2pcf-omp", stem+"3pcf-omp"] if geometry == "1d" and statistic == "both"
                     else [stem+suffix+"-omp"])
            detail = ("Preserve integer forest IDs and the intended same-forest exclusions. The 3D LOS-tree hierarchy can use certified cell/range moments "
                      "(--lya2-kernel 1, --lya3-kernel 5) while retaining zero bin slops. Measure performance on your forest lengths and separation range. "
                      "A 1D radial statistic and a 3D angular/LOS statistic are different estimators. Finite Legendre reconstruction also requires a multipole-convergence check.")
        elif geometry == "3d":
            names = ["octree-3pcf-3d-omp"]
            detail = ("This is the physical 3D Legendre-multipole estimator, distinct from angular kappa/shear statistics. "
                      "For survey counts use the supported data/random estimator and its random catalogue, weights and geometry. "
                      "Confirm the required 2PCF/3PCF outputs in the capability entry before running.")
        else:
            return Answer(question, "Please identify the estimator as well as the catalogue: angular scalar correlations, spin-2 shear, physical 3D scalar/data-random multipoles, or Lyα forest correlations.",
                          ["docs/search_methods.rst"], context=dict(self.context))
        available = self.available()
        verified = [name for name in names if name in available]
        text = detail+"\n\n"
        if verified:
            text += "Compatible candidates in this executable:\n"+"\n".join("  "+name for name in verified)
        else:
            text += "Source-documented candidates (not verified in this executable):\n"+"\n".join("  "+name for name in names)
        text += ("\n\nBenchmark elapsed time at the same measured approximation penalty. "
                 "Compare raw, normalized and corrected observables against an exact reference with the same data, bins and multipoles. "
                 "A subsample error does not certify the full catalogue. /topic measured-penalty explains the comparison commands.")
        if "perlmutter" in str(self.context.get("machine", "")).lower():
            text += ("\n\nOn Perlmutter, first measure one node, then ranks x threads with explicit core binding. "
                     "MPI replicates catalogues/trees and adds reductions; it is not automatically faster. "
                     "Keep approximation controls identical on all ranks. The CPU engines do not automatically use GPUs. "
                     "Independent realizations can often run as separate scheduler jobs. /topic mpi explains the limits.")
        return Answer(question, text, ["docs/search_methods.rst", "capabilities/engines.json", "docs/APPROXIMATION_DEFAULTS.md"],
                      topic="engine-selection", context=dict(self.context), notice=self.build_notice)

    def topic(self, name):
        found = next((p for p in self.knowledge.passages if p.kind == "answer" and p.topic == name), None)
        if found is None:
            raise ValueError("Unknown topic. Use /topics for the available names.")
        return self.from_passage(name, found)

    def from_passage(self, question, passage):
        notice = ""
        text = passage.text
        sources = self.knowledge.references(passage)
        if passage.stale:
            notice = "Some source files changed since this saved explanation was reviewed. Check the current excerpts below before using its settings or claims."
            current = self.knowledge.search(question, limit=3, kind="document")
            text += "\n\nCurrent documentation excerpts:\n"+"\n\n".join(p.text[:1500] for _, _, p in current)
            sources += [p.citation for _, _, p in current]
        return Answer(question, text, sources, topic=passage.topic, context=dict(self.context), notice=notice)

    def ask(self, question):
        if not question.strip() or len(question) > 12000:
            raise ValueError("Please use a question between 1 and 12,000 characters.")
        self.infer_context(question)
        lowered = question.lower().strip()
        query = question
        followup = bool(re.match(r"^(and |what about|how about|does it|is it|why is that|can it|then |those |that |same )", lowered))
        if followup and self.last:
            fresh = self.knowledge.search(question, limit=1, kind="answer")
            if not fresh or fresh[0][0] < 5 or fresh[0][1] < .5:
                query = question+" "+self.last.question
        changing_field = (followup and self.last and self.last.topic == "engine-selection"
                          and set(terms(question)) & {"kappa", "shear", "lya", "counts"})
        if changing_field or re.search(r"\b(which|what)\s+(?:(?:searching|shear|kappa|lya|spherical|angular|3d|1d|omp|mpi)\s+){0,3}(engine|method)\b|\b(best|appropriate|recommend|choose|select).{0,70}(engine|method)|engine.{0,30}(best|appropriate|recommend)", lowered):
            answer = self.recommend(question)
        elif lowered in {"hello", "hi", "hola", "help"}:
            answer = Answer(question, "Ask about an engine, approximation controls, normalization, MPI, timing or plots. /topics lists the saved scientific explanations; /help lists commands.")
        else:
            matches = self.knowledge.search(query, limit=6)
            saved = self.knowledge.search(query, limit=2, kind="answer")
            if saved and saved[0][0] >= 5 and saved[0][1] >= .30:
                answer = self.from_passage(question, saved[0][2])
            elif matches and matches[0][0] >= 4 and matches[0][1] >= .3:
                selected = [p for _, _, p in matches[:3]]
                answer = Answer(question, "I found these relevant passages; I do not have a reviewed answer for this precise question.\n\n"+
                    "\n\n".join(p.title+"\n"+p.text[:2200] for p in selected),
                    [p.citation for p in selected], context=dict(self.context))
            else:
                answer = Answer(question, "I could not find enough local evidence to answer reliably. Please name the engine, statistic or parameter, or use /search with its exact name. /topics shows the questions already covered.",
                                context=dict(self.context))
            if self.model and matches and any(not p.stale for _, _, p in matches):
                # Fresh documents are supplied before saved answers when the
                # latter were reviewed against a different source revision.
                selected = [p for _, _, p in matches if not p.stale][:5]
                try:
                    answer.text = self.model.answer(question, selected, self.context, self.history)
                    answer.sources = [p.citation for p in selected]
                    answer.mode = "local-model"
                    answer.notice = "Local-model synthesis; citations below identify the supplied evidence, not an automatic verification of every generated claim."
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    answer.notice = "Local model unavailable; showing offline help. "+str(exc)[:250]
        self.last = answer
        self.history.append(dict(question=question, answer=answer.text, sources=answer.sources, mode=answer.mode))
        return answer
