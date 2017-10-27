maxThreadFileBufferSize:
There are 3 conditions:

- maxThreadFileBufferSize == 0 or on all processors the local data size
is larger than the buffer:
No threading used; all files written immediately. Master will block
until whole file is written (so all processors have sent their data).
No limit on file size (since processors only send their data once the
master processor needs to write their chunk)

- total size of data (on all processors) <= maxThreadFileBufferSize:
Master waits until there is space in the thread buffer. Data is collected
on master and stored in the  write thread buffer.

- total size of data (on all processors) > maxThreadFileBufferSize
Master waits until there is space (for the local data only!) in the thread
buffer. Data gets stored in the write thread buffer. The write thread
picks up the request and collects the data processor by processor when
needed for writing (so there is no limit on the file size)
