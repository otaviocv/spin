import pytest as pt

from ..spin import SPIN


class TestSPINArguments:

    not_allowed_args = [
            pt.param("asdfjhsadf", marks=pt.mark.xfail),
            pt.param("sfs", marks=pt.mark.xfail),
            pt.param("srs", marks=pt.mark.xfail),
            pt.param("sys", marks=pt.mark.xfail),
            pt.param("sgs", marks=pt.mark.xfail),
            pt.param("neigborhood", marks=pt.mark.xfail),
            pt.param("neighborhoods", marks=pt.mark.xfail)
            ]

    @pt.mark.parametrize("argument", not_allowed_args)
    def test_error_raise_in_method_argument(self, argument):
        SPIN(method=argument)

    allowed_args = [
            pt.param("sts"),
            pt.param("neighborhood"),
            ]

    @pt.mark.parametrize("argument", allowed_args)
    def test_correct_string_in_method_argument(self, argument):
        SPIN(method=argument)
