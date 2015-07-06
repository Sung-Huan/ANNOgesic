import unittest
import os
import sys
import shutil
sys.path.append(".")
from annogesiclib.controller import Controller
from annogesiclib.projectcreator import ProjectCreator
from annogesiclib.paths import Paths

class ArgMock(object):

    subcommand = "test"
    project_path = "a_test_project"
    test_folder = os.path.join(project_path, "test_folder")

    def setup_folders(self, subcommand, test_folder):
        return {subcommand: [test_folder] + \
                self.required_test_folders(test_folder)}

    def required_test_folders(self, test_folder):
        return [os.path.join(test_folder, "test1"), 
                os.path.join(test_folder, "test2")]


class TestController(unittest.TestCase):

    def setUp(self):
        arg_mock = ArgMock()
        self.subcommand = arg_mock.subcommand
        self.test_folder = arg_mock.test_folder
        self.controller = Controller(arg_mock)
        self.test_project_name = arg_mock.project_path
        self.setup_folders = arg_mock.setup_folders(
                             self.subcommand, self.test_folder)
        self._version = 0.1
        self.create_subfolders = ProjectCreator().create_subfolders

    def tearDown(self):
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)

    def test_subcommand(self):
        self.controller.create_project(self._version)
        self.create_subfolders(self.setup_folders[self.subcommand])
        self.assertEqual(
        set(list(os.listdir(os.path.join(self.test_folder)))),
        set(['test1', 'test2']))

    def test_create_project(self):
        self.controller.create_project(self._version)
        self.assertEqual(
        set(list(os.listdir(self.test_project_name))),
        set(['input', 'output', 'used_annogesic_version.txt']))

if __name__ == "__main__":
    unittest.main()
