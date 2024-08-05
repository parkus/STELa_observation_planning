import warnings

import numpy as np


def midpts(ary, axis=None):
    """Computes the midpoints between points in a vector.

    Output has length len(vec)-1.
    """
    if type(ary) != np.ndarray: ary = np.array(ary)
    if axis == None:
        return (ary[1:] + ary[:-1])/2.0
    else:
        hi = np.split(ary, [1], axis=axis)[1]
        lo = np.split(ary, [-1], axis=axis)[0]
        return (hi+lo)/2.0


def mids2edges(mids, start='mid', first='adjacent', simple=False):
    """
    Reconstructs bin edges given only the midpoints.

    Parameters
    ----------
    mids : 1-D array-like
        A 1-D array or list of the midpoints from which bin edges are to be
        inferred.
    start : {'left'|'right'|'mid'}, optional
        left : start by assuming the spacing between the first two midpts
            is the same as the spacing between the first two bin edges and
            work from there
        right : same as above, but using the last two midpoints and working
            backwords
        mid : put one bin edge at the middle of the center two midpoints and
            work outwards. If there is an odd number of midpoints, use the
            middle midpoint and the one to the left to start.
    first : {float|'adjcacent'|'linear-i'|'linear-x'|function}, optional
        Width of the starting bin to the spacing between the midpoints to extrapolate the
        width of the first or last bin.

        The linear options try to extrapolate the width of the start bin by
        assuming the bin widths follow the same linear change as the midpoint
        spacings. linear-i' assumes a linear change with respect to the bin index
        whereas 'linear-x' assumes a linear change with respect to the bin
        value. These can only be used with start set to
        'left' or 'right'. Note that using 'linear' can produce nonsensical output
        if the spacing between midpoints does not vary linearly.

        Alternatively, a function may be provided that fits the midpoint spacings.
        That function should
        return a fit function (just like scipy's interp1d does), such that if
        result = function(xin,yin), then yout = result(xout) is a valid.

    Result
    ------
    edges : np.array
        The inferred bin edges.

    Could be accelerated with a cython implementation.
    """

    if simple:
        edges = midpts(mids)
        d0 = edges[0] - mids[0]
        d1 = mids[-1] - edges[-1]
        return np.insert(edges, [0, len(edges)], [mids[0] - d0, mids[-1] + d1])

    mids = np.array(mids)
    N = len(mids)
    e = np.zeros(N+1)
    if type(first) is not float and first != 'adjacent' and start == 'mid':
        raise ValueError("Start can only be 'mid' if fit == 'none'.")

    if type(first) is float:
        if start == 'left': e[0] = mids[0] - first/2.0
        if start == 'right': e[-1] = mids[-1] + first/2.0
    elif first == 'adjacent':
        if start == 'left': e[0] = mids[0] - (mids[1] - mids[0])/2.0
        if start == 'right': e[-1] = mids[-1] + (mids[-1] - mids[-2])/2.0
    else:
        d = mids[1:] - mids[:-1]
        x = midpts(mids)
        if first == 'linear-x':
            c = np.polyfit(x, d, 1)
            fitfun = lambda x: np.polyval(c, x)
        if first == 'linear-i':
            cdi = np.polyfit(np.arange(N-1), d, 1)
            cix = np.polyfit(x, np.arange(N-1), 2)
            def fitfun(x):
                i = np.polyval(cix, x)
                return np.polyval(cdi, i)
        elif callable(first):
            fitfun = first(x,d)

        if start == 'left':
            d0 = fitfun(mids[0])
            e[0] = mids[0] - d0/2.0
        if start == 'right':
            d1 = fitfun(mids[-1])
            e[-1] = mids[-1] + d1/2.0

    if start == 'left':
        for i in np.arange(0,N): e[i+1] = 2*mids[i] - e[i]
    if start == 'right':
        for i in np.arange(N-1,-1,-1): e[i] = 2*mids[i] - e[i+1]
    if start == 'mid':
        i = N//2
        e[i] = (mids[i-1] + mids[i])/2.0
        for i in np.arange(i-1,-1,-1): e[i] = 2*mids[i] - e[i+1]
        for i in np.arange(i+1,N): e[i+1] = 2*mids[i] - e[i]

    if any(e[1:] - e[:-1] <= 0):
        warnings.warn('There are zero or negative length bins in the output. '
                      'Consider using a different fit or start.', RuntimeWarning)

    return e


def cumtrapz(y, x, zero_start=False):
    result = np.cumsum(midpts(y)*np.diff(x))
    if zero_start:
        result = np.insert(result, 0, 0)
    return result


def intergolate(x_bin_edges,xin,yin, left=None, right=None):
    """Compute average of xin,yin within supplied bins.

    This funtion is similar to interpolation, but averages the curve repesented
    by xin,yin over the supplied bins to produce the output, yout.

    This is particularly useful, for example, for a spectrum of narrow emission
    incident on a detector with broad pixels. Each pixel averages out or
    "dilutes" the lines that fall within its range. However, simply
    interpolating at the pixel midpoints is a mistake as these points will
    often land between lines and predict no flux in a pixel where narrow but
    strong lines will actually produce significant flux.

    left and right have the same definition as in np.interp
    """

    x = np.hstack((x_bin_edges, xin))
    x = np.sort(x)
    y = np.interp(x, xin, yin, left, right)
    I = cumtrapz(y, x, True)
    Iedges = np.interp(x_bin_edges, x, I)
    y_bin_avg = np.diff(Iedges)/np.diff(x_bin_edges)
    return y_bin_avg



def common_axes(fig, subplotspec=111, covering_axs=None):
    if covering_axs is None:
        bigax = fig.add_subplot(subplotspec)
    else:
        positions = [ax.get_position() for ax in covering_axs]
        x0 = min(p.x0 for p in positions)
        y0 = min(p.y0 for p in positions)
        x1 = max(p.x1 for p in positions)
        y1 = max(p.y1 for p in positions)
        pos = [x0,y0,x1-x0,y1-y0]
        bigax = fig.add_axes(pos)
    [bigax.spines[s].set_visible(False) for s in ['top', 'bottom', 'left', 'right']]
    bigax.tick_params(labelleft=False, labelbottom=False)
    bigax.tick_params(which='both', left=False, bottom=False)
    bigax.set_zorder(-10)
    return bigax


def safe_interp_table(xnew, oldxcol, oldycol, table):
    usable = np.ones(len(table), bool)
    for col in [oldxcol, oldycol]:
        if hasattr(table[col], 'mask'):
            usable = usable & ~table[col].mask
    slim = table[usable]
    slim.sort(oldxcol)
    return np.interp(xnew, slim[oldxcol], slim[oldycol])