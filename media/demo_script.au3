#include <Constants.au3>

;
; AutoIt Version: 3.0
; Language:       English
; Platform:       Win9x/NT
; Author:         Jonathan Bennett (jon at autoitscript dot com)
;
; Script Function:
;   Opens Notepad, types in some text and then quits the application.
;

; Prompt the user to run the script - use a Yes/No prompt with the flag parameter set at 4 (see the help file for more details)
Local $iAnswer = MsgBox(BitOR($MB_YESNO, $MB_SYSTEMMODAL), "AutoIt Example", "This script will run Notepad, type in some text and then quit.  Do you want to run it?")

; Check the user's answer to the prompt (see the help file for MsgBox return values)
; If "No" was clicked (7) then exit the script
If $iAnswer = 7 Then
	MsgBox($MB_SYSTEMMODAL, "AutoIt", "OK.  Bye!")
	Exit
EndIf

; Run Notepad
Run("powershell.exe")

; Wait for the Notepad to become active. The classname "Notepad" is monitored instead of the window title
WinWaitActive("Windows PowerShell")
sleep(500)

; Now that the Notepad window is active type some text
Send("cd C:\Users\Parza\Dropbox\University\Code\finite_field_arithmetic {ENTER}")
sleep(2000)

Opt('SendKeyDelay', 100);
Send("python {ENTER}")
Sleep(500)
Send("from finite_field_arith import * {ENTER}")
Sleep(500)
Send("F2 = Field(2, 1, None){ENTER}")
Send("irr = Polynomial([1, 1, 0, 1], F2){ENTER}")
Send("print(irr){ENTER}")
Send("F8 = Field(2, 3, irr){ENTER}")
Send("print(F8){ENTER}")
Send("f = Polynomial([1, 0, 1], F2){ENTER}")
Send("f_in_field = Element(f, F8){ENTER}")
Send("print(f_in_field){ENTER}")
Send("generated = f_in_field.generated_subgroup(){ENTER}")
Send("for elem in generated:{ENTER}")
Send("    print(elem){ENTER}{ENTER}")

; Now wait for Notepad to close before continuing
WinWaitClose("Windows PowerShell")

; Finished!
