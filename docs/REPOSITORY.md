# Repository maintenance policy

The active implementation lives in `collaborator_rebuild`. The previous
release remains available in Git history and through `docs/released_study.md`.
Moving exploratory artifacts does not retract or validate their claims.

Track source code, analysis settings, provenance manifests, concise status
documents, and selected evidence needed to reproduce a claimed result.
Keep downloaded papers, private feedback, environments, caches, raw archives,
and regenerable figure collections local. Never commit collaborator messages
or third party publications merely because they exist in the workspace.

Current derived active outputs are ignored because they remain exploratory.
Promote a compact, documented evidence set only after the scientific settings
and provenance are resolved. Existing committed release evidence is preserved.

The active instruction DOCX and thesis PDF remain local root files and are
ignored. Their content is summarized in the analysis specification, rather
than republished. Downloads needed for the support audit are obtained by
`python -m collaborator_rebuild.check_support` with provenance recorded locally.

Remaining root scripts and numbered paper figures belong to the previous
release or its reproduction tools. They remain in place because their paths
are referenced by manuscript and pipeline workflows. Moving those paths
requires a separate dependency migration and release verification.

No history rewrite, remote deletion, or push is part of this organization.
The reorganization archives obsolete exploratory files rather than destroying
potentially useful evidence. Ignored local archives are left untouched.
