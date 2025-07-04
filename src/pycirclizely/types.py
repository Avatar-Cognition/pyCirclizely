from typing import (
    Any,
    Callable,
    List,
    Literal,
    Optional,
    Sequence,
    Union,
)

import numpy as np

Numeric = Union[int, float]
NumericSequence = Union[Sequence[Numeric], np.ndarray]
NumericComponent = Union[int, float, NumericSequence]

HoverText = Optional[Union[List[str], Literal["default"]]]

LabelFormatter = Optional[Callable[[float], str]]
TextFormatter = Optional[Callable[[str], str]]
HoverTextFormatter = Optional[Callable[[Any], List[str]]]
