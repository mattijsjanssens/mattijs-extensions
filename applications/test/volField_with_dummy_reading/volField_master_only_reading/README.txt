Problems
--------
- src/lagrangian/basic/Cloud/CloudIO.C:bool Foam::Cloud<ParticleType>::writeObject
Cloud only gets output if there is size.
Detection (headerOk etc) only gets done if there is a lagrangian/<cloudName>
on the local processor.
- src/meshTools/searchableSurface/triSurfaceMesh.C:bool Foam::triSurfaceMesh::writeObject
Uses OFstream directly.
- src/parallel/distributed/distributedTriSurfaceMesh/distributedTriSurfaceMesh.C:bool Foam::distributedTriSurfaceMesh::writeObject
Uses OFstream directly.

- uniform/time gets written to processorXXX, not processors


- regIOobject bit:
    - headerOk(syncPar, localValid) with optional typename checking
    - filePath(syncPar, localValid)
    - readHeaderOk
    - read
    - writeObject? But gets objectPath - maybe writeFilePath?

- masterOFstream and masterCollectingOFstream could be combined by
passing in a 'writer' that does something with the strings.

- switch between scheduled and nonblocking inside master*FileOperation?

- type in header is decomposedBlockData
- lives in processors/ subdirectory
- writing gets done to list of fileHandlers
- reading (of objects, not generic files) gets done from a list
  of fileHandlers. Start-up is
    (localFileHandler)
  and then from e.g. argList becomes
    (masterCollatingFileHandler masterFileHandler)

- system/controlDict gets read without

- needs to be done in thread

- is there a need for slaves to read their own slice?

- running distributed data?

- handle timeStampMaster. How do we know low-level down not to use
parallel version.

- headerOk() + readData() in one function so we don't open file twice

- handle local-only files? Threading? Store request, gather them and
  execute them. Or store and wait until next synchronised call is
  and we gather them then. Might be overkill to have thread polling all the
  time.
