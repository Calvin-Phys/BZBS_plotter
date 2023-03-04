import sys
import mechanize
import ssl
from bs4 import BeautifulSoup

def cif2json(fileName):
    """
    This function converts .cif file to .json file using the online tool 'SeeK-Path'
    """

    # open seekpath website
    ssl._create_default_https_context = ssl._create_unverified_context
    br = mechanize.Browser()
    br.set_handle_robots(False)
    br.set_handle_equiv(False)
    br.open('https://tools.materialscloud.org/seekpath/')
#     response = br.response()
    
    # fill form and submit
    br.select_form(enctype = 'multipart/form-data')
    br.form.add_file(open(fileName),content_type="application/octet-stream",name="structurefile")
    br.form["fileformat"] = ["cif-ase"]
    response2 = br.submit()

    # prase the response
    soup = BeautifulSoup(response2.read(), "html.parser")
    # print(soup.title)
    
    json_data = soup.find("code", {"id": "rawcodedata"})
    soup2 = BeautifulSoup(json_data.get_text(), "html.parser")
    
    f = open(fileName.replace('.cif','.json'), "w")
    f.write("".join(soup2.get_text().split()))
    f.close()

if __name__ == "__main__":
    cif2json(sys.argv[1])



