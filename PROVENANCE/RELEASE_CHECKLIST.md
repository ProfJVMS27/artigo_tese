# Public-release checklist

Do not create a public GitHub release or Zenodo archive until every applicable item is completed.

- [ ] Confirm that all coauthors approve the exact archived contents and version.
- [ ] Choose and add an explicit code license and a data license.  No license is implied by this package.
- [ ] Confirm that no solver license, proprietary executable, institutional credential, personal contact detail, or restricted operational information is included.
- [ ] Review every public input deck or model parameter file separately before adding it.  This package intentionally omits solver binaries and transient solver output.
- [ ] Create the GitHub repository from this directory, preserving the directory structure and `PROVENANCE/sha256_manifest.csv`.
- [ ] Tag the submission-ready release in GitHub.
- [ ] Connect the GitHub repository to Zenodo, create the archived release, and record the concept DOI and version DOI in the manuscript and repository citation file.
- [ ] Update author names, affiliations, ORCIDs, funding, and the data/code availability statement in the manuscript.
- [ ] Independently clone the repository into a new directory and verify the hashes and documented analysis entry points.
