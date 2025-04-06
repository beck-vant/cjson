import numpy as np
import json


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.int64):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        return super().default(obj)


def dump(mol: dict):
    json_dump = json.dumps(mol, sort_keys=False, indent=2, cls=NumpyEncoder)
    return json_dump


def load(cjson_str: str):
    return json.loads(cjson_str)
