from conans import ConanFile, CMake

class HelloConan(ConanFile):
    name = "libfixkalman"
    version = "20161008"
    url = "https://github.com/sunsided/libfixkalman.git"
    license = "MIT"
    author = "Markus Mayer (widemeadows@gmail.com)"
    requires = "libfixmath/20141230@sunside/stable", "libfixmatrix/20140117@sunside/stable"
    generators = "cmake"
    settings = "os", "compiler", "build_type", "arch"
    exports = "*.c", "*.h", "CMakeLists.txt", "*.rst", "AUTHORS", "LICENSE"

    def build(self):
        cmake = CMake(self.settings)
        self.run('cmake %s %s' % (self.conanfile_directory, cmake.command_line))
        self.run("cmake --build . %s" % cmake.build_config)

    def package(self):
        self.copy("*.h", dst="include")
        self.copy("*.c", dst="src")
        self.copy("*.lib", dst="lib", src="lib")
        self.copy("*.a", dst="lib", src="lib")

    def package_info(self):
        name = "libfixkalman"
        self.cpp_info.libs = ["libfixkalman"]