"""A changed numerical acceptance definition must invalidate both published tables."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import generate_accuracy_docs as docs


def test_accuracy_documents_match_cases_and_detect_changed_tolerance(monkeypatch):
    expected=docs.generated()
    for path,content in expected.items():
        assert path.read_text()==content
    cases=[dict(c) for c in docs.SCALAR]
    cases[0]['rtol']=.0031
    monkeypatch.setattr(docs,'SCALAR',cases)
    for path,content in docs.generated().items():
        assert '0.0031' in content
        assert content!=expected[path]
