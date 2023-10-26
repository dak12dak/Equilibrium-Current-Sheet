@echo off 
IF EXIST "main.exe" ( IF EXIST "bg.exe" ( del "bg.exe" ) )
IF EXIST "main.exe" ( IF EXIST "bin\bg.exe" ( del "bin\bg.exe" ) )
IF EXIST "main.exe" ( REN "main.exe" "bg.exe" & move "bg.exe" "bin\bg.exe"  )
IF EXIST "main.o"   ( DEL "main.o" )
