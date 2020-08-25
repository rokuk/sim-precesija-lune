import utils


# Čas v dnevih
start = 0
stop = 365250

# Št delilnih intervalov
N = (stop-start)*24

# Reši in shrani
utils.savestate(utils.solve(start, stop, N))
