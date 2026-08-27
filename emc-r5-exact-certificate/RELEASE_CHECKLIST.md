# Release checklist

- [ ] Obtain agreement from all authors on the public license and replace the
      provisional `LICENSE` file.
- [ ] Create the GitHub repository and add its canonical URL to
      `CITATION.cff` and `README.md`.
- [ ] Run both verification commands in a fresh Python 3.12.13 environment.
- [ ] Confirm that `git status --short` remains clean after verification.
- [ ] Regenerate `SHA256SUMS.txt` after every release-file change.
- [ ] Create an immutable release tag such as `v1.0.0`.
- [ ] Build a certificate-only release archive and record its SHA-256 digest.
- [ ] If an archival DOI is minted, add it to `CITATION.cff` and the paper's
      data-availability or reproducibility statement.
- [ ] Verify that the paper cites the immutable release rather than only the
      mutable default branch.
