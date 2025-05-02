import numpy as np
from typing import ClassVar, Tuple, List
from pycirclizely import config

class PolarSVGPatchBuilder:
    """Builds Plotly-compatible SVG paths by approximating arcs with line segments."""
    
    points: ClassVar[float] = config.ARC_POINTS

    @staticmethod
    def _polar_to_cart(theta: float, r: float) -> Tuple[float, float]:
        """Convert polar to Cartesian coordinates with Plotly orientation."""
        adjusted_theta = -(theta - np.pi / 2)  # 0°=up, clockwise
        x = r * np.cos(adjusted_theta)
        y = r * np.sin(adjusted_theta)
        return (x, y)
    
    @classmethod
    def arc_rectangle(cls, radr: Tuple[float, float], width: float, height: float) -> str:
        """Create rectangular arc sector approximated with line segments."""
        min_rad, min_r = radr
        max_rad = min_rad + width
        max_r = min_r + height

        # Bottom arc points (clockwise)
        bottom_points = [
            cls._polar_to_cart(angle, min_r)
            for angle in np.linspace(min_rad, max_rad, cls.points)
        ]

        # Top arc points (counter-clockwise)
        top_points = [
            cls._polar_to_cart(angle, max_r)
            for angle in np.linspace(max_rad, min_rad, cls.points)
        ]

        # Build path
        path = f"M {bottom_points[0][0]} {bottom_points[0][1]}"
        for point in bottom_points[1:]:
            path += f" L {point[0]} {point[1]}"
        for point in top_points:
            path += f" L {point[0]} {point[1]}"
        path += " Z"
        
        return path

    @staticmethod
    def _svg_path_from_points(
        points: list[tuple[float, float]], closed: bool = False
    ) -> str:
        if not points:
            return ""
        path = [f"M {points[0][0]} {points[0][1]}"]
        path += [f"L {x} {y}" for x, y in points[1:]]
        if closed:
            path.append("Z")
        return " ".join(path)

    @classmethod
    def arc_line(cls, rad_lim: Tuple[float, float], r_lim: Tuple[float, float]) -> str:
        """Create smooth arc between two points."""
        rad_start, rad_end = rad_lim
        r_start, r_end = r_lim

        if rad_start == rad_end:
            # Straight radial line
            start = cls._polar_to_cart(rad_start, r_start)
            end = cls._polar_to_cart(rad_end, r_end)
            return f"M {start[0]} {start[1]} L {end[0]} {end[1]}"

        # Generate points along the arc
        angles = np.linspace(rad_start, rad_end, cls.points)
        radii = np.linspace(r_start, r_end, cls.points)
        points = [cls._polar_to_cart(angle, radius) 
                 for angle, radius in zip(angles, radii)]

        # Build path
        path = f"M {points[0][0]} {points[0][1]}"
        for point in points[1:]:
            path += f" L {point[0]} {point[1]}"
            
        return path

    @classmethod
    def line(cls, rad_lim: Tuple[float, float], r_lim: Tuple[float, float]) -> str:
        """Create straight line between two points."""
        start = cls._polar_to_cart(rad_lim[0], r_lim[0])
        end = cls._polar_to_cart(rad_lim[1], r_lim[1])
        return f"M {start[0]} {start[1]} L {end[0]} {end[1]}"

    @classmethod
    def arc_arrow(cls, rad: float, r: float, drad: float, dr: float,
                 head_length: float = np.pi/90, shaft_ratio: float = 0.5) -> str:
        """Create directional arrow with arc segments."""
        shaft_size = dr * shaft_ratio
        y_shaft_bottom = r + ((dr - shaft_size) / 2)
        y_shaft_upper = r + dr - ((dr - shaft_size) / 2)
        
        is_forward = drad >= 0
        drad = abs(drad)
        head_length = min(head_length, drad)
        
        rad_shaft_tip = rad + (drad - head_length) if is_forward else rad - (drad - head_length)
        rad_arrow_tip = rad + drad if is_forward else rad - drad
        
        # Calculate key points
        p1 = cls._polar_to_cart(rad, y_shaft_bottom)
        p2 = cls._polar_to_cart(rad_shaft_tip, y_shaft_bottom)
        p3 = cls._polar_to_cart(rad_shaft_tip, r)
        p4 = cls._polar_to_cart(rad_arrow_tip, (r + r + dr)/2)
        p5 = cls._polar_to_cart(rad_shaft_tip, r + dr)
        p6 = cls._polar_to_cart(rad_shaft_tip, y_shaft_upper)
        p7 = cls._polar_to_cart(rad, y_shaft_upper)
        
        # Build path
        path = f"M {p1[0]} {p1[1]}"
        
        # Bottom shaft (arc segment)
        angles = np.linspace(rad, rad_shaft_tip, cls.points//3)
        for angle in angles[1:]:
            point = cls._polar_to_cart(angle, y_shaft_bottom)
            path += f" L {point[0]} {point[1]}"
        
        # Arrow head
        path += f" L {p3[0]} {p3[1]}"
        path += f" L {p4[0]} {p4[1]}"
        path += f" L {p5[0]} {p5[1]}"
        
        # Top shaft (arc segment)
        angles = np.linspace(rad_shaft_tip, rad, cls.points//3)
        for angle in angles[1:]:
            point = cls._polar_to_cart(angle, y_shaft_upper)
            path += f" L {point[0]} {point[1]}"
        
        path += " Z"
        return path