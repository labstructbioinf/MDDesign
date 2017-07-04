def prError(prt): print("\033[91m{}\033[00m" .format(prt))
def prOK(prt): print("\033[92m{}\033[00m" .format(prt))
def prWarning(prt): print("\033[93m{}\033[00m" .format(prt))
def prHeader(prt): print("\033[94m {}\033[00m" .format(prt))
def get_fn(model, rep, step=0, type=None):
    fn = "%s/%s_%s/" % (model, model, rep)
    if type == 'traj':
        assert(step != 0)
        fn += "_%s/" % step
        fn += '%s_%s_%s.nc' % (model, rep, step)
    elif type == 'traj_ai':
        assert(step != 0)
        fn += "_%s/" % step
        fn += '%s_%s_%s_ai.nc' % (model, rep, step)
    elif type == 'rst':
        if step == 0:
            fn += "%s_%s_i_s.rst7" % (model, rep)
        else:
            fn += "_%s/" % step
            fn += '%s_%s_%s.rst7' % (model, rep, step)
    elif type == 'mdg':
        fn += "%s_%s.mdg" % (model, rep)
    elif type == 'top':
        fn += "%s_%s.parm7" % (model, rep)
    elif type == 'solv_top':
        fn += "%s_%s_i_s.parm7" % (model, rep)
    elif type == 'inp':
        assert (step != 0)
        fn += "_%s/" % step
        fn += '%s_%s_%s.in' % (model, rep, step)
    elif type == 'sh':
        assert (step != 0)
        fn += "_SH/%s.sh" % step
    return fn
