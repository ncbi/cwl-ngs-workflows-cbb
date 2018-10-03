import unittest
import cwltool


class TestCWL(unittest.TestCase):

    def test_cwl(self):
        fac = cwltool.factory.Factory()
        echo = fac.make("echo.cwl")
        result = echo(inp="foo")
        print(result)


if __name__ == '__main__':
    unittest.main()
