import pytest

from CRISPRSCope.editing_rate_ci import (
    _depth_stability_layout,
    _interval_axis_limits,
    _row_plot_height,
)


def test_interval_plot_height_scales_with_amplicon_count():
    assert _row_plot_height(1) == pytest.approx(3.2)
    assert _row_plot_height(1, legend=True) == pytest.approx(3.8)
    assert _row_plot_height(30) > _row_plot_height(5)
    assert _row_plot_height(30, legend=True) > _row_plot_height(1, legend=True)


def test_interval_axis_limits_focus_on_plotted_confidence_intervals():
    editing_limits = _interval_axis_limits(
        [97.7, 98.2],
        [97.9, 98.7],
        minimum_span=5.0,
        bounds=(0.0, 100.0),
    )
    assert editing_limits == pytest.approx((95.0, 100.0))

    effect_limits = _interval_axis_limits(
        [0.4],
        [0.9],
        minimum_span=2.0,
        include_zero=True,
    )
    assert effect_limits[0] < 0
    assert effect_limits[1] > 0.9
    assert effect_limits[1] - effect_limits[0] == pytest.approx(2.0)


@pytest.mark.parametrize(
    "n_amplicons, expected_rows, expected_columns, expected_width",
    [
        (1, 1, 1, 8.0),
        (2, 1, 2, 12.0),
        (3, 1, 3, 18.0),
        (4, 2, 2, 12.0),
        (5, 2, 3, 18.0),
        (30, 10, 3, 18.0),
    ],
)
def test_depth_stability_layout_scales_to_amplicon_count(
    n_amplicons,
    expected_rows,
    expected_columns,
    expected_width,
):
    rows, columns, width, height = _depth_stability_layout(n_amplicons)
    assert (rows, columns) == (expected_rows, expected_columns)
    assert width == pytest.approx(expected_width)
    assert height == pytest.approx(1.5 + (3.8 * expected_rows))


def test_depth_stability_layout_rejects_empty_plot():
    with pytest.raises(ValueError, match="At least one amplicon"):
        _depth_stability_layout(0)
