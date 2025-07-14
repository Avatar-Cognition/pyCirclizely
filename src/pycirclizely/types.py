from typing import (
    Any,
    Callable,
    Literal,
    Optional,
    Sequence,
    TypeAlias,
    Union,
)

import numpy as np

Numeric: TypeAlias = Union[int, float]
NumericSequence: TypeAlias = Union[Sequence[Numeric], np.ndarray]
NumericComponent: TypeAlias = Union[int, float, NumericSequence]

HoverText: TypeAlias = Optional[Union[list[str], Literal["default"]]]

LabelFormatter: TypeAlias = Optional[Callable[[float], str]]
TextFormatter: TypeAlias = Optional[Callable[[str], str]]
HoverTextFormatter: TypeAlias = Optional[Callable[[Any], list[str]]]
