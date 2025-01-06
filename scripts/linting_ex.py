'''Example script to demonstrate linting'''

# Inputs
inventory = "GFAS"
ver = "v1.2"
region = "global"
res = "0.1"
species = 'co'
year = 2024

load_file_path = "C:/Users/vlqc7666/OneDrive - University of Leeds/Documents/PhD/data/GFAS" #load files from here
save_file_path = "C:/Users/vlqc7666/OneDrive - University of Leeds/Documents/PhD/project_wetland/data/GFAS" #path where you want to save the files

filename = f"GFAS_v1.2_{species}_{year}.nc"
load_file = f"{load_file_path}/{filename}"
print(f"{year} is standardised.")
months = [f"{month:02d}" for month in range(1, 11)] #change depending on end date

for month_counter, month in enumerate(months):
    save_file = f"{save_file_path}/{year}/{inventory}_{ver}_{region}_{res}_{species}_{year}{month}.nc"

    # loading and standardising the data
    print(f"{month}/{year} is loaded.")

    # saving the file
    print(f"{month}/{year} is saved.")
