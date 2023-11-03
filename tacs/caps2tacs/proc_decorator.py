__all__ = ["root_proc", "parallel", "root_broadcast"]

from functools import wraps


# do something on all active procs
def parallel(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        if self.comm.rank in self.active_procs:
            return method(self, *args, **kwargs)
        else:

            def empty_function(self):
                return self  # for potential method cascading

            return empty_function

    return wrapped_method


# Define a root proc decorator so that certain ESP/CAPS method like pre and post analysis
# are only ran on the root processor
def root_proc(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        if self.comm.rank == self.active_procs[0]:
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
            root = self.active_procs[0]
            if self.comm.rank == root:
                output = method(self, *args, **kwargs)
            output = self.comm.bcast(output, root=root)
            return output

    return wrapped_method
