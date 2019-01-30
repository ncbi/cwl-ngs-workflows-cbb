import unittest
import cwltool
import cwltool.factory


class TestSet(unittest.TestCase):

    def test_cwl(self):
        fac = cwltool.factory.Factory()

        echo = fac.make("./tools/basic/echo.cwl")
        result = echo(stdout="echo.stdout",msg="foo")
        self.assertEqual(result['out_stdout']['basename'], 'echo.stdout')
