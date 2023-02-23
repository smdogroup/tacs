__all__ = ["root_proc"]

from functools import wraps


# Define a root proc decorator so that certain ESP/CAPS method like pre and post analysis
# are only ran on the root processor
def root_proc(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        if self.comm is None or self.comm.rank == 0:
            return method(self, *args, **kwargs)
        else:

            def empty_function(self):
                return

            return empty_function

    return wrapped_method


# define a decorator to broadcast certain outputs to other processors
def root_broadcast(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        if self.comm is None:
            return method(self, *args, **kwargs)
        else:
            output = None
            if self.comm.rank == 0:
                output = method(self, *args, **kwargs)
            output = self.comm.bcast(output, root=0)
            return output

    return wrapped_method
