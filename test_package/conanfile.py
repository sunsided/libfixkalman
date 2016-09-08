from conans import ConanFile, CMake
import os

channel  = os.getenv("CONAN_CHANNEL", "stable")
username = os.getenv("CONAN_USERNAME", "sunside")

class HelloReuseConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = "libfixkalman/20161008@%s/%s" % (username, channel)
    generators = "cmake"

    def build(self):
        cmake = CMake(self.settings)
        self.run('cmake "%s" %s' % (self.conanfile_directory, cmake.command_line))
        self.run("cmake --build . %s" % cmake.build_config)

    def test(self):
        self.run(os.sep.join([".","bin", "package-test"]))