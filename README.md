# CATNIP-Protoplanetary-Centering-and-Deprojection
This is python code that will take a fits image of a disk and find the center and inclination returning a deprojected image of the disk. 
README: 
Hello! Welcome to the current iteration of the deprojection and centering code for the CATNIP project! 
The goal of this work is to have a computer find both the true center and the deprojection for a protoplanetary disk fits image. 
Currently the inspiration code from Kate Follette is kept in the IDL file center_cirlesym.pro 
The new code that is written in python is found in protoplan_center.py 

To get started using this there is a fits file uploaded in the images folder of disk SGR V4046. However, I would recommend checking if that is still up to date. 
To run this code on your computer you will need a working conda. 

To Run: 
First, clone this repository and navigate to the folder in terminal:

cd Protoplanet Center + Deprojection
Then you will need to update the images' file path. In protoplan_center.py, on line 413, there is a line folder = Path('path/images'). Please replace this with your file path. 
Then create and activate a new virtual environment:

python -m venv venv
source venv/bin/activate # Windows: venv\Scripts\activate

Install the required dependencies:

pip install -r requirements.txt

and run the code: 

python protoplan_center.py 

This should ideally show you the center and deprojection for disk SGR V4046. If for some reason it doesn't: by all means try another fits image from ALMA, SPHERE, or Subaru, 
or reach out to me on github. 
