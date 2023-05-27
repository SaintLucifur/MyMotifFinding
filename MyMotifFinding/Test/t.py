header = """
<html> <head> <title> Output Data in an HTML file
</title> </head> <h1>MMF Known Motif Enrichment Results {dir}</h1>
<h2>TOP 10 Motif Found for <u>dataset</u></h2>
<style>
</style>
<body>
<table style="width:100%">
    <tr>
        <th>Motif</th>
        <th>Scores</th>
""".format(dir="ME")
def main():
    print(header.format(dir="Charles"))
    
main()