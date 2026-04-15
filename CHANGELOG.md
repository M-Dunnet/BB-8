# CHANGE LOG
## BB-8 (fork of C3POa) - Modified version
## 2026-04-16 - Current working version

### Added
- Added splint UMI detection, which allows BB-8 to determine chimeric reads. Reads with multiple UMIs across its subreads will be split into multiple seperate reads.
- BB-8 post-processing now adds read directionality to read headers so we know which strand was seqeunced. Ideal for strand bias in variant calling.

## Changed
- generateConsensus was returned to a function, rather than its own script that is called via subprocess. This requires more RAM but runs faster.
- abpoa executable dependency was removed; returned to pyabpoa. This means less writing temp files to disk; requires more RAM but runs faster
- racon executable dependency was removed; smoothing was slowing down processing and did not appear to do much. 
- Removed live consensus calling; ideally will add this back later but not needed for current applications
- removed 'no splint' options; not needed for current applications
- general code tidy; added comments to describe whats happening.
