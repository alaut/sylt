import matplotlib.pyplot as plt

from sylt.geometry import bivariate_binomial
from sylt.plotting import plot_phase_space

for mu in [0.5, 1, 1.12, 1.5, 2, 5]:

    x, y = bivariate_binomial(2, 1, 1_000_000, mu=mu)

    plot_phase_space(
        {'x': x, 'y': y},
        [['x', 'y']],
        title=f"{mu}",
    )

plt.show()

print('done!')
