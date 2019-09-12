from .sprim import SPrim

class SArm(SPrim):
    """
    Scattering kernel for Arm component - does not alter the neutron state.

    Parameters
    ----------
    None

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, idx=0, ctx=0, **kwargs):
        return

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        return