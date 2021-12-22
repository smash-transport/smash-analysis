import numpy as np

def lorentz_boost(four_vector, three_velocity):
    """Apply a Lorentz boost to a given 4-vector.

    Returns a tuple with the boosted 4-vector. This is designed so it works
    nicely with numpy."""
    assert len(four_vector) == 4
    assert len(three_velocity) == 3
    x = four_vector
    u = three_velocity
    if 0 == u[0] == u[1] == u[2]:
        return x
    uu = u[0]**2 + u[1]**2 + u[2]**2
    gamma = 1 / np.sqrt(1 - uu)
    k = (gamma - 1) / uu
    return (
        gamma * (x[0] - u[0]*x[1] - u[1]*x[2] - u[2]*x[3]),
        -gamma*u[0]*x[0] + (k*u[0]**2 + 1)*x[1] + k*(u[0]*u[1]*x[2] + u[0]*u[2]*x[3]),
        -gamma*u[1]*x[0] + (k*u[1]**2 + 1)*x[2] + k*(u[0]*u[1]*x[1] + u[1]*u[2]*x[3]),
        -gamma*u[2]*x[0] + (k*u[2]**2 + 1)*x[3] + k*(u[0]*u[2]*x[1] + u[1]*u[2]*x[2]),
    )


def test_lorentz_boost():
    r = [10., 1., 0., 0.]
    u = [0.5, 0, 0]
    r_boosted = lorentz_boost(r, u)
    r_expected = [10.9697, -4.6188, 0., 0.]
    assert all(abs(a - b) < 6e-5 for a, b in zip(r_boosted, r_expected))


if __name__ == '__main__':
    test_lorentz_boost()
    print('All tests finished successfully.')
