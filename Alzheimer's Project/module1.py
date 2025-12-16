#opening and returning header/row information for both datasets
import csv 
with open("UpdatedMetaData.csv", newline='') as f: 
    reader = csv.reader(f) 
    headers = next(reader) 
    print('The demographic headers are:') 
    for h in headers: 
        print(h) 
    for r in reader:
        pass
    print('The totals rows is: %d' % reader.line_num) 

with open("UpdatedLuminex.csv", newline="") as f2:
    reader = csv.reader(f2) 
    headers = next(reader) 
    print('The Luminex headers are:') 
    for h in headers: 
        print(h) 
    for r in reader:
        pass
    print('The totals rows is: %d' % reader.line_num) 