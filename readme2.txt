open terminal--->run the script as "python.lscode1.py"//the script name
once you run this the user interface asks to put the path of the log file
once we input the path, the script asks 

Choose what to analyze:
1) lp-bp* interactions
2) bp-bp* interactions
Enter your choice (1, 2 ):

so choose any one between 1 or 2 
after we input the number it asks 

if we choose 1 it further asks to choose the below type of interactions 
Choose interaction type:
a) n-σ* interactions
b) n-π* interactions
Enter a or b:

if we choose 2 then it shows 
Choose interaction type:
c) σ-σ* interactions
d) σ-π* interactions
e) π-π* interactions
f) π-σ* interactions
Enter c, d, e or f:

so choose the appropriate type of interactions u need to extract 

it asks user to set the threshold energy ,once you enter the threshold energy the output file will have those particular type of interaction  above the threshold level set by the user. 
Enter minimum energy threshold:

after choosing the type of interactions the results will be saved in the final_output2.txt which contains only the particular type of interaction information tht the user has choosen.

i ran the code for N2 conformer by taking NO2_N2.log file as input ,after running the script i got 4 output files. They were interactions2.txt,extracted_data2.txt,lncN2.txt and final_output2.txt.

interactions2.txt contains information about all the type of interactions .

extracted_data2.txt contains information about the extracted data from log file that is  bond hybridization data and second order pertubation table 

lncN2.txt contains the information about the type of interaction along with the second order pertubation table data and the shortest path between the donor and the first acceptor atom.

final_output2.txt  contains only the particular type of interaction having the energy which are more than the threshold energy set by user.


