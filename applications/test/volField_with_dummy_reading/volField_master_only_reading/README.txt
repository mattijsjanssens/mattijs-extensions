- decomposedBlockData

- system/controlDict gets read without

- is fileHandler not server.

- needs to be done in thread

- needs to use scheduled communication optionally so master
  does not need to read all of file.

- change file format to be List<List<char>> so we can read without
  index file.

- is there a need for slaves to read their own slice?

- running distributed data?

- handle timeStampMaster. How do we know low-level down not to use
parallel version.

- headerOk() + readData() in one function so we don't open file twice

- findInstance master-only version

- handle local-only files? Threading? Store request, gather them and
  execute them. Or store and wait until next synchronised call is
  and we gather them then. Might be overkill to have thread polling all the
  time.

- processor-indexed regIOobjects.
    readStream(const word& typeName)
    - master looks for objectPath()
    - checks type:
        - volScalarField : read as is
        - dictionary :
            - check for objectPath() "+.index" file
                (- read size and check if agrees with)
                - lookup 'processorXXX' to give start. Problem: start is
                  64 bits. Encode as hex?
                - fseek to start
                - check that starts with key?
            - else
                - subDict 'processorXXX'
            - headerOk() returns headerOk from processor0?
            - or type is Dictionary<volScalarField> and typeName is
              volScalarField.

    - regIOobject::writeObject
        - are we using indexed writing?
            no: use server().ofstream() etc.
        - yes: send to master (server().ofstream())
        - Open 'dictionary' file
        - insert contributions as subdicts with 'processorXXX' as key
        - write objectPath + ".index" with starting points of subdicts
