"""Scientific tutor behavior: no survey calculation or model installation needed."""
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import socket
import sqlite3
import subprocess
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tutor.catalog import inspect_catalog
from tutor.cli import command, run_report, save_session
from tutor.knowledge import Knowledge
from tutor.local_model import LocalModel, local_endpoint, NoRedirect
from tutor.session import Tutor


@pytest.fixture(scope="module")
def knowledge():
    k = Knowledge(ROOT, ":memory:")
    yield k
    k.close()


@pytest.fixture
def tutor(knowledge, monkeypatch):
    # Offline retrieval must not depend on a reachable service.
    def forbidden(*args, **kwargs):
        raise AssertionError("Offline tutor attempted a network connection")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    return Tutor(knowledge)


@pytest.mark.parametrize("question,topic", [
    ("How does cTreeBalls work?", "how-it-works"),
    ("How do I require 10 percent error for every nonzero coefficient?", "measured-penalty"),
    ("How do I profile scalar two-ball computation?", "profiling"),
    ("Why does --selection sofia-fig1 fail for 8 to 200 arcmin edges?", "binning"),
    ("How can I reduce memory pressure?", "catalog-memory"),
    ("How do I transfer catalog ownership?", "catalog-memory"),
    ("What is a pivot?", "how-it-works"),
    ("What are the default approximation parameters?", "defaults"),
    ("Does CBALLS_SHEAR_BIN_THETA=0.1 overwrite --tree-theta 1?", "separate-controls"),
    ("What do you mean by zero tolerance?", "zero-tolerance"),
    ("Will 0.5 give at most 5% relative L2 error?", "measured-penalty"),
    ("What is the role of CBALLS_SHEAR_PROFILE=1?", "profiling"),
    ("Does NMultipoles play a role in kappa and shear engines?", "multipoles"),
    ("Does normalization need all N_n count multipoles?", "normalization"),
    ("Can MPI make improvements?", "mpi"),
    ("What about running MPI on Perlmutter?", "perlmutter"),
    ("Which units does the spherical engine use?", "geometry"),
    ("What does resume do?", "resume"),
    ("Why is octree shear slower than balltree and kdtree?", "octree-speed"),
    ("Why is CPU time much greater than real time?", "cpu-time"),
    ("Why do plots ignore min-sep and max-sep?", "binning"),
])
def test_scientific_question_retrieval(tutor, question, topic):
    answer = tutor.ask(question)
    assert answer.topic == topic
    assert answer.sources and answer.mode == "offline"
    assert not answer.notice


def test_followups_and_engine_context(tutor):
    tutor.ask("What is CBALLS_SHEAR_PROFILE?")
    assert tutor.ask("Is it the same as verbose?").topic == "profiling"
    assert tutor.ask("And MPI?").topic == "mpi"
    answer = tutor.ask("Which engine is best for DES shear 3PCF?")
    assert answer.topic == "engine-selection"
    assert "balltree-shear-sphere" in answer.text
    assert "phase 3 and bin allowance 0.5" in answer.text
    assert "nsmooth defaults to 8" in answer.text
    answer = tutor.ask("And kappa?")
    assert answer.context["field"] == "kappa"
    assert "Scalar reuse" in answer.text
    assert "balltree-2balls-omp" in answer.text
    command(tutor, "/reset")
    assert not tutor.context and not tutor.history


def test_unknown_question_and_missing_context(tutor):
    answer = tutor.ask("What is the capital of France?")
    assert not answer.sources and "could not find enough" in answer.text
    answer = tutor.ask("Which engine should I use?")
    assert "specify field, geometry" in answer.text


def test_context_rejects_invalid_update_atomically(tutor):
    tutor.set_context('field=shear geometry=sphere catalog="DES sample"')
    before = dict(tutor.context)
    with pytest.raises(ValueError):
        tutor.set_context("field=kappa threads=-4")
    assert tutor.context == before


def test_references_have_current_sources(knowledge):
    answers = [p for p in knowledge.passages if p.kind == "answer"]
    assert len(answers) >= 29
    for passage in answers:
        assert not passage.stale, passage.topic
        assert passage.references
        for source in passage.references:
            assert source in knowledge.file_hashes, (passage.topic, source)
        assert passage.topic in (ROOT/passage.source).read_text().splitlines()[passage.line-1]


@pytest.fixture
def tiny_root(tmp_path):
    root = tmp_path/"source with spaces"
    (root/"docs").mkdir(parents=True)
    (root/"tutor/knowledge").mkdir(parents=True)
    doc = root/"docs/sample.md"
    doc.write_text("# Spherical geometry\n\nUse chord distances for sphere bins.\n")
    question = dict(id="geometry", title="Spherical geometry", questions=["How are sphere bins defined?"],
        keywords=["chord"], answer="Use chord distances.", sources=["docs/sample.md"],
        source_sha256={"docs/sample.md": hashlib.sha256(doc.read_bytes()).hexdigest()})
    (root/"tutor/knowledge/questions.json").write_text(json.dumps({"questions": [question]}, indent=2))
    return root


def test_index_refresh_add_remove_and_stale_answer(tiny_root):
    k = Knowledge(tiny_root)
    assert k.updated
    k.close()
    k = Knowledge(tiny_root)
    assert not k.updated
    k.close()
    (tiny_root/"docs/sample.md").write_text("# Spherical geometry\n\nUpdated chord definition.\n")
    (tiny_root/"docs/new.md").write_text("# Forests\nKeep forest identifiers.\n")
    k = Knowledge(tiny_root)
    assert k.updated
    assert any(p.stale for p in k.passages)
    assert k.search("forest identifiers")[0][2].source == "docs/new.md"
    k.close()
    (tiny_root/"docs/new.md").unlink()
    k = Knowledge(tiny_root)
    assert "docs/new.md" not in k.file_hashes
    k.close()


def test_symlinks_do_not_import_outside_files(tiny_root, tmp_path):
    private = tmp_path/"private.md"
    private.write_text("private test data")
    (tiny_root/"docs/link.md").symlink_to(private)
    k = Knowledge(tiny_root)
    assert "docs/link.md" not in k.file_hashes
    k.close()


@pytest.mark.parametrize("cache_kind", ["directory", "corrupt", "readonly-refresh"])
def test_unusable_cache_falls_back(tiny_root, tmp_path, cache_kind):
    cache = tmp_path/"database"
    if cache_kind == "directory":
        cache.mkdir()
    elif cache_kind == "corrupt":
        cache.write_text("not a sqlite database")
    else:
        k = Knowledge(tiny_root, cache)
        k.close()
        cache.chmod(0o444)
        (tiny_root/"docs/sample.md").write_text("Changed geometry information.")
    try:
        k = Knowledge(tiny_root, cache)
        assert k.database == ":memory:"
        assert k.passages and k.cache_notice
        k.close()
    finally:
        if cache.is_file():
            cache.chmod(0o644)


def test_stale_answer_includes_current_passages(knowledge, monkeypatch):
    passage = next(p for p in knowledge.passages if p.topic == "separate-controls")
    monkeypatch.setattr(passage, "stale", True)
    answer = Tutor(knowledge).topic("separate-controls")
    assert "changed" in answer.notice
    assert "Current documentation excerpts" in answer.text


def test_compiled_availability_and_unverified_fallback(tutor, tmp_path):
    if not (ROOT/"cballs").exists():
        pytest.skip("Build cballs for native availability checks")
    available = tutor.available()
    assert available and "Availability checked" in tutor.build_notice
    tutor.executable = tmp_path/"absent"
    tutor._available = None
    assert not tutor.available()
    tutor.set_context("field=shear geometry=sphere")
    assert "not verified" in tutor.recommend().text


def test_recommendation_respects_geometry_and_lya_1d(tutor):
    tutor.set_context("field=shear geometry=flat")
    assert "Do not pass flat" in tutor.recommend().text
    tutor.set_context("field=lya geometry=sphere")
    answer = tutor.recommend()
    assert "insufficient" in answer.text
    assert all((ROOT/p).is_file() for p in answer.sources)
    tutor.set_context("geometry=1d statistic=both")
    answer = tutor.recommend()
    assert "lya-1d-tree-2pcf" in answer.text
    assert "lya-1d-tree-3pcf" in answer.text


def npz_headers(path, columns):
    with zipfile.ZipFile(path, "w") as archive:
        for name, shape, dtype in columns:
            header = repr(dict(descr=dtype, fortran_order=False, shape=shape)).encode()+b"\n"
            # Deliberately omit array payload: successful inspection proves it
            # does not deserialize object data or read full survey columns.
            archive.writestr(name+".npy", b"\x93NUMPY\x01\x00"+len(header).to_bytes(2,"little")+header)


def test_npz_headers_only_and_context_reset(tutor, tmp_path):
    path = tmp_path/"catalog with spaces.npz"
    npz_headers(path, [("gamma1", (5_000_000,), "<f8"), ("gamma2", (5_000_000,), "<f8"),
                      ("metadata", (), "|O"), ("ra", (5_000_000,), "<f8"), ("dec", (5_000_000,), "<f8")])
    data = inspect_catalog(path)
    assert data["rows"] == 5_000_000 and data["fields"] == ["shear"]
    assert data["geometry_hint"] == "sphere"
    command(tutor, f'/catalog "{path}"')
    assert tutor.context["field"] == "shear"
    npz_headers(path, [("unrecognized", (17,), "<f8")])
    command(tutor, f'/catalog "{path}"')
    assert "field" not in tutor.context and "geometry" not in tutor.context


def test_bad_npz_does_not_end_conversation(tutor, tmp_path):
    path = tmp_path/"bad.npz"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("bad.npy", b"\x93NUMPY")
    with pytest.raises(ValueError, match="Invalid NPY"):
        inspect_catalog(path)


def test_fits_headers_optional(tmp_path):
    fits = pytest.importorskip("astropy.io.fits")
    path = tmp_path/"catalog.fits"
    table = fits.BinTableHDU.from_columns([fits.Column(name=n, format="D", array=[1.,2.])
                                          for n in ("RA", "DEC", "KAPPA")])
    table.writeto(path)
    data = inspect_catalog(path)
    assert data["rows"] == 2 and data["fields"] == ["kappa"]
    assert data["geometry_hint"] == "sphere"


def test_saved_penalty_and_unmeasured_are_distinct(tutor, tmp_path):
    path = tmp_path/"approximation-penalty.json"
    path.write_text(json.dumps(dict(status="UNMEASURED", target_percent=5, observables={})))
    assert "does not mean zero error" in run_report(tutor, path).text
    data = dict(status="MEASURED", target_percent=5, native_wall_speedup=1.6,
                observables={"shear": dict(relative_l2_percent=6, relative_target_met=False, absolute_l2=.02),
                             "zero_norm": dict(relative_l2_percent=None, relative_target_met=False, absolute_l2=0)})
    path.write_text(json.dumps(data))
    text = run_report(tutor, tmp_path).text
    assert "6%; above target" in text and "undefined" in text
    assert "1.6x" in text and "does not establish" in text


def test_transcript_explicit_and_no_overwrite(tutor, tmp_path):
    tutor.ask("What is zero tolerance?")
    assert not list(tmp_path.iterdir())
    path = tmp_path/"session.json"
    save_session(tutor, path)
    data = json.loads(path.read_text())
    assert data["turns"][0]["sources"] and data["turns"][0]["mode"] == "offline"
    with pytest.raises(FileExistsError):
        save_session(tutor, path)


@pytest.mark.parametrize("url", ["https://example.org", "http://127.0.0.1.evil.example", "ftp://localhost", "http://localhost/path", "http://user:pass@localhost", "http://localhost?x=1", "http://localhost:invalid"])
def test_model_endpoint_rejects_nonlocal_or_invalid(url):
    with pytest.raises(ValueError):
        local_endpoint(url)


@pytest.mark.parametrize("url", ["http://localhost:11434", "http://127.0.0.1:11434/", "http://[::1]:11434"])
def test_model_endpoint_accepts_loopback(url):
    assert local_endpoint(url).endswith("/api/chat")


def test_model_protocol_context_sources_and_history(tutor, monkeypatch):
    model = LocalModel("already-installed-test-model")
    captured = {}
    def open_response(request, timeout):
        captured.update(json.loads(request.data))
        assert request.full_url == "http://127.0.0.1:11434/api/chat"
        return io.BytesIO(json.dumps(dict(message=dict(content="Measured errors need a reference [1]."))).encode())
    monkeypatch.setattr(model.opener, "open", open_response)
    tutor.ask("What is zero tolerance?")
    tutor.model = model
    answer = tutor.ask("What are the default approximation parameters?")
    assert answer.mode == "local-model" and answer.sources
    assert captured["stream"] is False and captured["model"] == "already-installed-test-model"
    messages = captured["messages"]
    assert "Retrieved source material:" in messages[0]["content"]
    assert messages[1]["content"] == "What is zero tolerance?"
    assert messages[-1]["role"] == "user"
    assert "not an automatic verification" in answer.notice


@pytest.mark.parametrize("failure", ["offline", "malformed", "empty", "oversized", "invalid-shape"])
def test_unavailable_model_falls_back_without_losing_answer(tutor, monkeypatch, failure):
    model = LocalModel("installed-test")
    def response(*args, **kwargs):
        if failure == "offline":
            raise OSError("connection refused")
        contents = {"malformed": b"not json", "empty": b'{"message":{"content":""}}',
                    "oversized": b"x"*2_000_001, "invalid-shape": b"[]"}
        return io.BytesIO(contents[failure])
    monkeypatch.setattr(model.opener, "open", response)
    tutor.model = model
    answer = tutor.ask("What is zero tolerance?")
    assert answer.mode == "offline" and answer.topic == "zero-tolerance"
    assert "Local model unavailable" in answer.notice


def test_no_redirects_or_environment_proxy(monkeypatch):
    monkeypatch.setenv("HTTP_PROXY", "http://example.org:8080")
    model = LocalModel("test")
    assert not any(type(h).__name__ == "ProxyHandler" and h.proxies for h in model.opener.handlers)
    with pytest.raises(ValueError, match="redirected"):
        NoRedirect().redirect_request(None, None, 302, "moved", {}, "http://example.org")


@pytest.fixture
def native_env():
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("CBALLS_TUTOR_") or key in {"OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MV2_COMM_WORLD_SIZE", "SLURM_PROCID", "SLURM_NTASKS"}:
            env.pop(key)
    env.update(CBALLS_TUTOR_PYTHON=sys.executable, PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1")
    return env


def native(args, env, cwd, stdin=None):
    binary = ROOT/"cballs"
    if not binary.is_file():
        pytest.skip("Build cballs to exercise the native launcher")
    return subprocess.run([str(binary), "options=IA-tutor", *args], cwd=cwd, env=env,
                          input=stdin, text=True, capture_output=True, timeout=30)


def test_native_launch_elsewhere_json_no_calculation(native_env, tmp_path):
    process = native(["--ask", "What is zero tolerance?", "--json", "--database", ":memory:"], native_env, tmp_path)
    assert process.returncode == 0, process.stderr
    answer = json.loads(process.stdout)
    assert answer["topic"] == "zero-tolerance"
    assert not list(tmp_path.iterdir()), "Tutor started writing calculation outputs"


def test_native_piped_conversation_and_eof(native_env, tmp_path):
    process = native(["--database", ":memory:"], native_env, tmp_path,
        "/context field=shear geometry=sphere statistic=3pcf\nWhat is CBALLS_SHEAR_PROFILE?\nIs it the same as verbose?\n/sources\n")
    assert process.returncode == 0, process.stderr
    assert "offline documentation mode" in process.stdout
    assert "CBALLS_SHEAR_PROFILE" in process.stdout and "Sources:" in process.stdout


def test_native_argument_forwarding_spaces_and_no_shell(native_env, tmp_path):
    root = tmp_path/"tutor resources with spaces"
    (root/"scripts").mkdir(parents=True)
    (root/"scripts/ia_tutor.py").write_text("import json,sys; print(json.dumps(sys.argv[1:]))\n")
    native_env["CBALLS_TUTOR_ROOT"] = str(root)
    question = 'literal $(touch unexpected-file) `touch unexpected-too` ; spaces'
    process = native(["--", "--ask", question], native_env, tmp_path)
    assert process.returncode == 0, process.stderr
    args = json.loads(process.stdout)
    assert args[-2:] == ["--ask", question] and args[:2] == ["--root", str(root)]
    assert not (tmp_path/"unexpected-file").exists()
    assert not (tmp_path/"unexpected-too").exists()


@pytest.mark.parametrize("setting,expected", [
    (("CBALLS_TUTOR_ROOT", "/nonexistent/cballs-tutor"), "Cannot find scripts"),
    (("CBALLS_TUTOR_PYTHON", "/nonexistent/python"), "Cannot start tutor Python"),
    (("OMPI_COMM_WORLD_SIZE", "2"), "serial conversation"),
    (("PMI_SIZE", "2"), "serial conversation"),
])
def test_native_failure_paths(native_env, tmp_path, setting, expected):
    native_env[setting[0]] = setting[1]
    process = native([], native_env, tmp_path)
    assert process.returncode != 0
    assert expected in process.stderr
    assert not list(tmp_path.iterdir())


def test_native_mixed_options_are_rejected(native_env, tmp_path):
    process = subprocess.run([str(ROOT/"cballs"), "options=IA-tutor,shear-pivot-reuse"],
         cwd=tmp_path, env=native_env, text=True, capture_output=True, timeout=10)
    assert process.returncode == 2 and "separate information mode" in process.stderr


def test_python_entry_stdlib_only(native_env, tmp_path):
    process = subprocess.run([sys.executable, "-S", "-B", str(ROOT/"scripts/ia_tutor.py"),
                              "--ask", "What is zero tolerance?", "--database", ":memory:", "--json"],
                             cwd=tmp_path, env=native_env, text=True, capture_output=True, timeout=30)
    assert process.returncode == 0, process.stderr
    assert json.loads(process.stdout)["mode"] == "offline"


def test_malformed_report_is_readable_error(tutor, tmp_path):
    path = tmp_path/"bad-report.json"
    path.write_text(json.dumps(dict(target_percent=5, observables=[])))
    with pytest.raises(ValueError, match="Invalid approximation-penalty"):
        run_report(tutor, path)


def test_native_question_can_contain_mode_token(native_env, tmp_path):
    process = native(["--ask", "options=IA-tutor", "--json", "--database", ":memory:"], native_env, tmp_path)
    assert process.returncode == 0, process.stderr
    assert json.loads(process.stdout)["question"] == "options=IA-tutor"


def test_native_slurm_multi_task_rejected(native_env, tmp_path):
    native_env.update(SLURM_PROCID="0", SLURM_NTASKS="2")
    process = native([], native_env, tmp_path)
    assert process.returncode == 2 and "multiple Slurm tasks" in process.stderr


def test_direct_python_invalid_local_model_falls_back(native_env, tmp_path):
    # A closed loopback port exercises an actual refused local connection,
    # with no model or external network dependency.
    process = native(["--ask", "What is zero tolerance?", "--json", "--database", ":memory:",
                      "--backend", "ollama", "--model", "test-not-installed", "--model-url", "http://127.0.0.1:1", "--timeout", "1"], native_env, tmp_path)
    assert process.returncode == 0, process.stderr
    answer = json.loads(process.stdout)
    assert answer["topic"] == "zero-tolerance" and answer["mode"] == "offline"
    assert "Local model unavailable" in answer["notice"]
