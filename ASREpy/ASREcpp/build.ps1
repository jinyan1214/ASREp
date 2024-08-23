Set-Location .\src\build\

cmake .. -G "Visual Studio 17 2022" -A x64

cmake --build . --config Release --parallel 6

cmake --install .

Set-Location ..\..