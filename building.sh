# if you want to debug
bazel build  //... --compilation_mode=dbg
# Then you just gdb to the file you want
gdbtui bazel-bin/bin/..

# if you want final build
bazel build  //...