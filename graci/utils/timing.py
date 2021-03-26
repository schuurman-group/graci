"""
Timers and timing information for subroutines.

Allows for nested timers. Constructed to ensure no overlap between timers,
so thereotically sum(timers) = total time
"""
import time

timer_list = dict()
active_stack = []

class timed:
    """A wrapper/descriptor for timing functions and methods.

    This takes advantage of the fact that __call__ is called for a
    wrapped function and __get__ is called for a wrapped class
    method. Functions are automatically named [module name].[func name]
    and methods are named [class name].[method name].
    """
    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        modname = self.func.__module__.split('.')[-1]
        return self._run_func(modname + '.' +
                              self.func.__name__, *args, **kwargs)

    def __get__(self, obj, type=None):
        def hooked(*args, **kwargs):
            return self._run_func(obj.__class__.__name__ + '.' +
                                  self.func.__name__, obj, *args, **kwargs)
        return hooked

    def _run_func(self, name, *args, **kwargs):
        start(name)
        result = self.func(*args, **kwargs)
        stop(name)
        return result


class Timer:
    """A Timer object for a process 'name'."""
    def __init__(self, name):
        self.name       = name
        self.calls      = 0
        self.wall_time  = 0.
        self.cpu_time   = 0.
        self.wall_start = 0.
        self.cpu_start  = 0.
        self.wall_kids  = 0.
        self.cpu_kids   = 0.
        self.running    = False

    def start_call(self):
        """Start the wall and CPU times and track latest call."""
        self.calls     += 1
        self.running    = True
        self.wall_start = time.time()
        self.cpu_start  = time.clock()
        self.wall_kids  = 0.
        self.cpu_kids   = 0.

    def stop_call(self, cumulative=False):
        """Add to wall_time and cpu_time. Returns change in wall and
        CPU times."""
        d_wall = time.time()  - self.wall_start
        d_cpu  = time.clock() - self.cpu_start

        # If cumulative, don't subtract off time spent in nested timers.
        if cumulative:
            self.wall_time += d_wall
            self.cpu_time  += d_cpu
        else:
            self.wall_time += d_wall - self.wall_kids
            self.cpu_time  += d_cpu - self.cpu_kids

        self.running = False
        return d_wall, d_cpu


def start(name):
    """Starts the timer 'name'.

    If 'name' doesn't exist, it is created and added to the list of timers.
    """
    global timer_list, active_stack

    # if timer already in list, activate it
    if name not in timer_list:
        timer_list[name] = Timer(name)

    # add the timer to the active stack
    active_stack.append(timer_list[name])

    if active_stack[-1].running:
        raise ValueError('START timer:' + str(name) +
                         ' called while running (could be recursion issue)')

    active_stack[-1].start_call()


def stop(name, cumulative=False):
    """Stops the timer 'name'.

    Trying to stop an unknown timer raises NameError. If not cumulative,
    the time of child functions are subtracted from the total run time.
    """
    global active_stack

    # for time being, assume, last function added to stack, is the first to
    # be stopped. If more flexibility is required, we'll address it at that time
    if active_stack[-1].name != name:
        raise NameError('STOP timer: ' + str(name) +
                        ' called, but timer not at top of active stack.\n')

    d_wall, d_cpu = active_stack[-1].stop_call(cumulative)

    # pop the child off the stack, then go through the active_stack and update
    # the child times for each of the parent timers
    active_stack.pop()

    # add the time for the child to the parent
    if len(active_stack) > 0:
        active_stack[-1].wall_kids += d_wall
        active_stack[-1].cpu_kids  += d_cpu


def print_timings():
    """Prints out a timing report, sorted from highest wall time to
    lowest.

    The timer 'global' is stopped automatically and its time is used as
    the cumulative total time.
    """
    global timer_list

    # ensure that the global timer has finished and get the total execution time
    if timer_list['global'].running:
        stop('global', cumulative=True)
    tot_cpu  = timer_list['global'].cpu_time
    tot_wall = timer_list['global'].wall_time

    sort_list = sorted(timer_list.items(),
                       key=lambda unsort: unsort[1].wall_time, reverse=True)

    # pass timing information as a string
    ostr =  '\n' + '-'*39 + ' timings summary ' + '-'*39 + ' \n'
    ostr += ('routine'.ljust(35) + 'calls'.rjust(12) +
             'wall time'.rjust(16) + 'frac.'.rjust(8) +
             'cpu time'.rjust(16)  + 'frac.'.rjust(8) + '\n')
    ofrm = '{:35s}{:12d}{:16.4f}{:8.2f}{:16.4f}{:8.2f}\n'
    frac_wall = 0.
    frac_cpu  = 0.
    for sort in sort_list:
        rout = str(sort[0])
        if rout != 'global':
            ncall = sort[1].calls
            wtim  = sort[1].wall_time
            ctim  = sort[1].cpu_time
            ostr += ofrm.format(rout, ncall, wtim, wtim/tot_wall, ctim,
                                ctim/tot_cpu)
            frac_wall += wtim/tot_wall
            frac_cpu  += ctim/tot_cpu

    ostr += '-'*95 + '\n'
    ostr += ('**total**'.ljust(47) +
             '{:16.4f}{:8.2f}{:16.4f}{:8.2f}\n\n'.format(tot_wall, frac_wall,
                                                             tot_cpu, frac_cpu))
    return ostr
