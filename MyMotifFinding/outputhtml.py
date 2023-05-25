def main():
    html = open("KnownMotifFinding.html", "w")
    header = """
<html>\n<head>\n<title> \nOutput Data in an HTML file \
</title>\n</head> <body><h1>Welcome to <u>GeeksforGeeks</u></h1>\
\n<h2>A <u>CS</u> Portal for Everyone</h2> \n</body></html>
"""
    dict = {"100":"Motif1", "101":"Motif2", "102":"Motif3"}
    html.write(header)
    for score in dict.keys():
        html.write("Motif: {0} has a score of {1}\n<br>".format(dict[score], score))
    html.close()

main()