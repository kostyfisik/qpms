import time
import sys

def _time_b(active = True, name = None, step = None):
    '''
    Auxiliary function for keeping track of elapsed time.
    Returns current time (to be used by _time_e).
    '''
    now = time.time()
    if active:
        if not name:
            name = sys._getframe(1).f_code.co_name
        if step:
            print('%.4f: %s in function %s started.' % (now, step, name), file = sys.stderr)
        else:
            print('%.4f: Function %s started.' % (now, name), file=sys.stderr)
        sys.stderr.flush()
    return now

def _time_e(start_time, active = True, name = None, step = None):
    now = time.time()
    if active:
        if not name:
            name = sys._getframe(1).f_code.co_name
        if step:
            print('%.4f: %s in function %s finished (elapsed %.2f s).' 
                    % (now, step, name, now - start_time), file = sys.stderr)
        else:
            print('%.4f: Function %s finished (elapsed %.2f s).' 
                    % (now, name, now - start_time), file = sys.stderr)
        sys.stderr.flush()

