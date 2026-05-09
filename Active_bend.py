import math

def make_coords(r=0.5, n=4, top_gap=1.8, angle_between_deg=30, y_top=0.0):
    """
    Generates center coordinates for two columns of circles arranged in a V-shape.
    """
    theta_deg = angle_between_deg / 2.0
    theta = math.radians(theta_deg)

    step = 2 * r
    ux = math.sin(theta)
    uy = -math.cos(theta)

    x_left0,  y_left0 = -top_gap / 2.0, y_top
    x_right0, y_right0 =  top_gap / 2.0, y_top

    left_col, right_col = [], []

    for i in range(n):
        xl = x_left0 - i * step * ux
        yl = y_left0 + i * step * uy
        xr = x_right0 + i * step * ux
        yr = y_right0 + i * step * uy
        left_col.append((xl, yl))
        right_col.append((xr, yr))

    return left_col[::-1] + right_col

def save_to_txt_simple(coords, filename="curve_150.txt"):
    """
    Saves coordinates in a 3-column tab-separated format: X  Y  Z
    """
    with open(filename, 'w') as f:
        for (x, y) in coords:
            # Format: X [tab] Y [tab] Z (fixed at 0.0)
            f.write(f"{x}\t{y}\t0.0\n")
    print(f"File successfully saved to: {filename}")

# --- Parameters ---
radius = 0.5
num_per_side = 3  # Total points = 150
gap = 1.7
angle = 30

# Generate and Save
coords_list = make_coords(r=radius, n=num_per_side, top_gap=gap, angle_between_deg=angle, y_top=0.0)
save_to_txt_simple(coords_list, "curve_150.txt")