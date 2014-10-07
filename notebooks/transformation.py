import numpy as np
from scipy import ndimage

def sum_inside_contours(index_field, index_levels, fields, area=1.):
    """Sum the arrays in args within countours
    of index_field. Bins are specifiied by index levels.
    index_levels must be monotonically increasing or decreasing."""
    di = np.diff(index_levels)
    sign = np.sign(di[0])
    assert np.all(np.sign(di)==sign)
    shape = index_field.shape
    labels = np.digitize(index_field.ravel(), index_levels,
                         right=(sign==1))
    labels.shape = shape
    
    res = []
    
    # wrap it in a list if we only got one argument
    if isinstance(fields, np.ndarray):
        fields = [fields,]
    for field in fields:
        assert field.shape == shape
        field = np.ma.masked_array(field, index_field.mask)
        res.append( ndimage.sum(
            (area*field).filled(0.),
            labels=labels, index=np.arange(len(index_levels))
        ))
    if len(res)==1:
        return res[0]
    else:
        return res