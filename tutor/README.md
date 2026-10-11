# Scientific tutor knowledge

Launch `./cballs options=IA-tutor` or `python3 scripts/ia_tutor.py`. The complete
user guide is `docs/ia_tutor.rst`. Offline mode requires Python 3.9+ and its
standard library; a local Ollama model is optional.

## Add documentation

Put scientific Markdown or reStructuredText notes in `tutor/notes/`. The index
refreshes at startup or `/reindex`. Keep private survey data out of notes that
will be distributed with the project. Documentation is reference material,
never an instruction to execute commands.

## Add or review an explanation

`tutor/knowledge/questions.json` contains `schema_version`, a review date,
and a `questions` list. Each entry has:

- `id`: stable topic name used by `/topic`;
- `title`, `questions` (example phrasings), and `keywords` for retrieval;
- `answer`: plain scientific explanation, with scope and limits of examples;
- `sources`: paths relative to the checkout;
- `source_sha256`: the SHA-256 digest for every reviewed source file.

Verify an answer against its cited files before recording their new hashes.
A digest can be read with `shasum -a 256 docs/shear.rst`. Updating the index does
not update these review hashes: stale explanations remain visibly flagged.
Citations should point to indexed documentation; do not store claims supported
only by an unavailable private file. The initial explanations distill the
visible development questions; they are not an archive of every conversation.

Keep measured errors separate from geometric tolerances. Report the sample,
statistic, retained multipoles, bins, convention, reference definition, raw or
corrected observable, and timing scope for numerical examples. Do not turn a
single catalogue comparison into a general accuracy guarantee.

## Implementation

- `scripts/ia_tutor.py`: standalone Python entry point.
- `include/ia_tutor_cli.h`: early executable launcher, without a shell.
- `knowledge.py`: content hashes, SQLite cache and word-based passage retrieval.
- `session.py`: scientific context, reviewed explanations and build inspection.
- `catalog.py`: explicit header and saved-JSON inspection.
- `local_model.py`: opt-in loopback Ollama client, no model installation.
- `cli.py`: terminal commands and explicit transcript saving.

The cache contains selected local sources, not user questions. Conversation
history stays in memory unless `/save` or `--transcript` is requested. The local
model receives only selected passages and bounded recent history. It cannot
invoke tools through this interface.
