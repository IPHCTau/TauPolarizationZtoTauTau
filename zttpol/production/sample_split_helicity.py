import law
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div
from columnflow.columnar_util import set_ak_column,remove_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)



@producer(
    uses={"helicity_sign"},
    produces={"process_id"},
    mc_only=True,
)
def split_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    # all process_ids from cmsdb ewk process, can be accessed automatically
    helpos_proc_id = 511236
    helneg_proc_id = 511235

    process_id = ak.where(events.helicity_sign < 0, helneg_proc_id, events.process_id)
    process_id = ak.where(events.helicity_sign > 0, helpos_proc_id, process_id)
    
    events = remove_ak_column(events, "process_id")
    events = set_ak_column(events, "process_id", process_id, value_type=np.int64)

    return events
