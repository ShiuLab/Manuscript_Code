@REM echo %time%
set "b=%CD%"
python %b%\image_converter.py %b%\images
@REM seed_images\template_testing
@chdir %b%\ij152-win-java8\ImageJ
@REM echo %time%
FOR %%a in (%b%\images\8bit_*) DO (
@REM echo %time%
java -cp ij.jar ij.ImageJ -macro %b%\small_plate_partial_macro.ijm %%a
@REM echo %time%
)
@REM echo %time%
@PAUSE