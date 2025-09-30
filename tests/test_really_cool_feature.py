import numpy as np
import pytest


@pytest.mark.parametrize("dtype", ["int8", "int16", "int32", "int64"])
def test_dtypes(dtype):
    assert str(np.dtype(dtype)) == dtype


@pytest.mark.parametrize(
    "dtype",
    [
        "float32",
        pytest.param("int16", marks=pytest.mark.skip),
        pytest.param("int32", marks=pytest.mark.xfail(reason="to show how it works")),
    ],
)
def test_mark(dtype):
    assert str(np.dtype(dtype)) == "float32"


@pytest.fixture
def fake_data():
    return np.array([1, 2, 3])


@pytest.fixture(params=["int8", "int16", "int32", "int64"])
def dtype(request):
    return request.param


def test_series(fake_data, dtype):
    result = fake_data.astype(dtype)
    assert result.dtype == dtype

    # expected = xr.DataArray(np.array([1, 2, 3], dtype=dtype))
    # assert_equal(result, expected)
