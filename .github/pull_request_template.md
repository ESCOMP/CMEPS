### Description of changes

### Specific notes

Contributors other than yourself, if any:

CMEPS Issues Fixed (include github issue #):

Are changes expected to change answers?
 - [ ] bit for bit
 - [ ] different at roundoff level
 - [ ] more substantial 

Any User Interface Changes (namelist or namelist defaults changes)?
 - [ ] Yes
 - [ ] No

Testing performed if application target is CESM:(either UFS-S2S or CESM testing is required):
- [ ] (required) CIME_DRIVER=nuopc scripts_regression_tests.py
   - machines:
   - details (e.g. failed tests):
- [ ] (required) CESM testlist_drv.xml
   - machines and compilers:
   - details (e.g. failed tests):
- [ ] (optional) CESM prealpha test
   - machines and compilers
   - details (e.g. failed tests):

Testing performed if application target is UFS-S2S:
- [ ] (required) UFS-S2S testing
   - description:
   - details (e.g. failed tests):

Hashes used for testing:
- [ ] CESM:
  - repository to check out: https://github.com/ESCOMP/CESM.git
  - branch: nuopc_dev
  - hash:
- [ ] UFS-S2S, then umbrella repostiory to check out and associated hash:
  - repository to check out:
  - branch:
  - hash:
