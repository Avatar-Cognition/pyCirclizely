from typing import ClassVar

import numpy as np

from pycirclizely import config


class PolarSVGPatchBuilder:
    """Builds Plotly-compatible SVG paths by approximating arcs with line segments."""

    n_points: ClassVar[int] = int(config.ARC_POINTS)  # Ensure this is an integer

    @staticmethod
    def _polar_to_cart(theta: float, r: float) -> tuple[float, float]:
        """Convert polar to Cartesian coordinates with Plotly orientation."""
        adjusted_theta = -(theta - np.pi / 2)  # 0°=up, clockwise
        x = r * np.cos(adjusted_theta)
        y = r * np.sin(adjusted_theta)
        return (x, y)

    @classmethod
    def arc_rectangle(
        cls, radr: tuple[float, float], width: float, height: float
    ) -> str:
        """Create rectangular arc sector approximated with line segments."""
        min_rad, min_r = radr
        max_rad = min_rad + width
        max_r = min_r + height

        # Bottom arc points (clockwise)
        bottom_points = [
            cls._polar_to_cart(float(angle), min_r)  # Explicit float conversion
            for angle in np.linspace(min_rad, max_rad, num=cls.n_points)
        ]

        # Top arc points (counter-clockwise)
        top_points = [
            cls._polar_to_cart(float(angle), max_r)  # Explicit float conversion
            for angle in np.linspace(max_rad, min_rad, num=cls.n_points)
        ]

        # Build path
        path = f"M {bottom_points[0][0]} {bottom_points[0][1]}"
        for point in bottom_points[1:]:
            path += f" L {point[0]} {point[1]}"
        for point in top_points:
            path += f" L {point[0]} {point[1]}"
        path += " Z"

        return path

    @classmethod
    def arc_line(cls, rad_lim: tuple[float, float], r_lim: tuple[float, float]) -> str:
        """Create smooth arc between two points."""
        rad_start, rad_end = rad_lim
        r_start, r_end = r_lim

        if rad_start == rad_end:
            # Straight radial line
            start = cls._polar_to_cart(rad_start, r_start)
            end = cls._polar_to_cart(rad_end, r_end)
            return f"M {start[0]} {start[1]} L {end[0]} {end[1]}"

        # Generate points along the arc
        angles = np.linspace(rad_start, rad_end, num=cls.n_points)
        radii = np.linspace(r_start, r_end, num=cls.n_points)
        points = [
            cls._polar_to_cart(float(angle), float(radius))  # Explicit float conversion
            for angle, radius in zip(angles, radii)
        ]

        # Build path
        path = f"M {points[0][0]} {points[0][1]}"
        for point in points[1:]:
            path += f" L {point[0]} {point[1]}"

        return path

    @classmethod
    def multi_segment_path(
        cls, rad: list[float], r: list[float], arc: bool = True
    ) -> str:
        """Builds a full line path by chaining segments between adjacent points."""
        if len(rad) != len(r):
            raise ValueError("rad and r must be the same length")

        segments = []
        for i in range(len(rad) - 1):
            segment = (
                cls.arc_line(rad_lim=(rad[i], rad[i + 1]), r_lim=(r[i], r[i + 1]))
                if arc
                else cls.straight_line(
                    rad_lim=(rad[i], rad[i + 1]), r_lim=(r[i], r[i + 1])
                )
            )
            # Remove the leading 'M' in all but the first segment, continuous path
            if i > 0 and segment.startswith("M"):
                segment = segment.replace("M", "L", 1)
            segments.append(segment)

        return " ".join(segments)

    @classmethod
    def straight_line(
        cls, rad_lim: tuple[float, float], r_lim: tuple[float, float]
    ) -> str:
        """Create straight line between two points."""
        start = cls._polar_to_cart(rad_lim[0], r_lim[0])
        end = cls._polar_to_cart(rad_lim[1], r_lim[1])
        return f"M {start[0]} {start[1]} L {end[0]} {end[1]}"

    @classmethod
    def arc_arrow(
        cls,
        rad: float,
        r: float,
        drad: float,
        dr: float,
        head_length: float = np.pi / 90,
        shaft_ratio: float = 0.5,
    ) -> str:
        """Create an arc-shaped arrow with proper path closure."""
        # Calculate shaft dimensions
        shaft_size = dr * shaft_ratio
        r_shaft_bottom = r + ((dr - shaft_size) / 2)
        r_shaft_upper = r + dr - ((dr - shaft_size) / 2)

        # Determine direction and adjust head length
        is_forward = drad >= 0
        abs_drad = abs(drad)
        head_length = min(head_length, abs_drad)

        # Calculate key positions
        if is_forward:
            rad_shaft_tip = rad + (abs_drad - head_length)
            rad_arrow_tip = rad + abs_drad
        else:
            rad_shaft_tip = rad - (abs_drad - head_length)
            rad_arrow_tip = rad - abs_drad

        # Build the bottom shaft
        bottom_shaft = cls.arc_line(
            rad_lim=(rad, rad_shaft_tip), r_lim=(r_shaft_bottom, r_shaft_bottom)
        )

        # Build the arrowhead
        arrow_bottom_tip = cls._polar_to_cart(rad_shaft_tip, r)
        arrowhead = cls._arrow_path(
            rad_start=rad_shaft_tip,
            angle=rad_arrow_tip,
            r_base=r + dr,
            r_tip=(r + (r + dr)) / 2,
        )
        arrow_path = f"L {arrow_bottom_tip[0]} {arrow_bottom_tip[1]} {arrowhead}"

        # Build the upper shaft (remove the 'M' command)
        upper_shaft = cls.arc_line(
            rad_lim=(rad_shaft_tip, rad), r_lim=(r_shaft_upper, r_shaft_upper)
        ).replace("M", "L", 1)

        path = f"{bottom_shaft} {arrow_path} {upper_shaft} Z"
        return path

    @classmethod
    def build_filled_path(
        cls,
        rad: list[float],
        r1: list[float],
        r2: list[float],
        arc: bool = True,
        closed: bool = True,
    ) -> str:
        """Build a closed path between two polar curves."""
        # Create upper path (from first to last point)
        upper_path = cls.multi_segment_path(rad, r1, arc=arc)

        # Create lower path (from last to first point)
        reversed_rad = list(reversed(rad))
        reversed_r2 = list(reversed(r2))
        lower_path = cls.multi_segment_path(reversed_rad, reversed_r2, arc=arc)

        # Combine paths (remove the initial 'M' from the lower path)
        if lower_path.startswith("M"):
            lower_path = "L" + lower_path[1:]

        # Combine and close if needed
        combined_path = f"{upper_path} {lower_path}"
        if closed:
            combined_path += " Z"

        return combined_path

    @classmethod
    def bezier_ribbon_path(
        cls,
        rad1_start: float,
        rad1_end: float,
        r1: float,
        rad2_start: float,
        rad2_end: float,
        r2: float,
        height_ratio: float,
        direction: int = 0,
        arrow_length_ratio: float = 0.05,
    ) -> str:
        """Create SVG path for a complete ribbon with proper arc connections."""
        # Calculate arrow radius adjustments
        arrow_r1 = r1 * (1 - arrow_length_ratio)
        arrow_r2 = r2 * (1 - arrow_length_ratio)

        # Handle different directions
        if direction in (-1, 2):
            # Start at arrow base on region1 and appply reverse arrow
            start_point = cls._polar_to_cart(rad1_start, arrow_r1)
            path = f"M {start_point[0]},{start_point[1]}"
            path += cls._arrow_path(rad1_end, (rad1_end + rad1_start) / 2, arrow_r1, r1)
        else:
            start_point = cls._polar_to_cart(rad1_start, r1)
            path = f"M {start_point[0]},{start_point[1]}"
            path += (
                " " + cls.arc_line((rad1_start, rad1_end), (r1, r1)).split(" ", 1)[1]
            )

        # Append Bezier curve to region2 (end1 → end2)
        path += " " + cls._bezier_segment(
            rad1_end,
            r1,
            rad2_end,
            arrow_r2 if direction in (1, 2) else r2,
            height_ratio,
        )

        # Handle forward arrow (pointing to middle of region2)
        if direction in (1, 2):
            path += cls._arrow_path(
                rad2_start, (rad2_start + rad2_end) / 2, arrow_r2, r2
            )
        else:
            path += (
                " " + cls.arc_line((rad2_end, rad2_start), (r2, r2)).split(" ", 1)[1]
            )

        # Append Bezier curve back to region1 (start2 → start1)
        path += " " + cls._bezier_segment(
            rad2_start,
            arrow_r2 if direction in (-1, 2) else r2,
            rad1_start,
            arrow_r1 if direction in (-1, 2) else r1,
            height_ratio,
        )

        # Close the path
        path += " Z"

        return path

    @classmethod
    def bezier_line_path(
        cls,
        rad1: float,
        r1: float,
        rad2: float,
        r2: float,
        height_ratio: float,
        direction: int = 0,
        arrow_height: float = 3.0,
        arrow_width: float = 0.05,
    ) -> str:
        """Create SVG path for a complete link line."""
        # Handle different directions
        if direction in (-1, 2):
            # Start at arrow base on region1 and appply reverse arrow
            arrow_r1 = r1 - arrow_height
            left_base = rad1 + arrow_width / 2
            right_base = rad1 - arrow_width / 2
            arrow_start = cls._polar_to_cart(left_base, arrow_r1)
            path = f"M {arrow_start[0]},{arrow_start[1]}"
            path += cls._arrow_path(
                right_base, (right_base + left_base) / 2, arrow_r1, r1
            )
            tip_point = cls._polar_to_cart(rad1, r1)
            path += f" M {tip_point[0]},{tip_point[1]}"

        else:
            start_point = cls._polar_to_cart(rad1, r1)
            path = f"M {start_point[0]},{start_point[1]}"

        # Append Bezier curve to region2 (end1 → end2)
        path += " " + cls._bezier_segment(rad1, r1, rad2, r2, height_ratio)

        # Handle forward arrow (pointing to middle of region2)
        if direction in (1, 2):
            arrow_r2 = r2 - arrow_height
            left_base = rad2 + arrow_width / 2
            right_base = rad2 - arrow_width / 2
            arrow_start = cls._polar_to_cart(left_base, arrow_r2)
            path += f" M {arrow_start[0]},{arrow_start[1]}"
            path += cls._arrow_path(
                right_base, (right_base + left_base) / 2, arrow_r2, r2
            )

        return path

    @classmethod
    def _bezier_segment(
        cls, rad1: float, r1: float, rad2: float, r2: float, height_ratio: float = 0.5
    ) -> str:
        """Generate SVG path for a Bezier curve segment between two points."""
        p2 = cls._polar_to_cart(rad2, r2)

        # Calculate control point
        mid_angle = (rad1 + rad2) / 2
        if height_ratio >= 0.5:
            r_ctl = (r1 + r2) * (height_ratio - 0.5)
            rad_ctl = mid_angle + np.pi  # Curve outward
        else:
            r_ctl = (r1 + r2) * (0.5 - height_ratio)
            rad_ctl = mid_angle  # Curve inward

        ctrl = cls._polar_to_cart(rad_ctl, r_ctl)

        return f" C {ctrl[0]},{ctrl[1]} {ctrl[0]},{ctrl[1]} {p2[0]},{p2[1]}"

    @classmethod
    def _arrow_path(
        cls, rad_start: float, angle: float, r_base: float, r_tip: float
    ) -> str:
        """Create properly oriented arrow at specified angle with proper closure."""
        # Calculate all points
        p_end = cls._polar_to_cart(rad_start, r_base)
        p_tip = cls._polar_to_cart(angle, r_tip)

        # Draw arrowhead and return to base point
        return (
            f" L {p_tip[0]},{p_tip[1]}"  # Draw to tip
            f" L {p_end[0]},{p_end[1]}"
        )  # Close arrow
