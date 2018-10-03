import unittest
import cwltool
import cwltool.factory


class TestSet(unittest.TestCase):

    def test_cwl(self):
        fac = cwltool.factory.Factory()

        echo = fac.make("./tools/echo.cwl")
        result = echo(inp="foo")
        print(result)
        self.assertEqual(result['out'], 'foo\n')
