import signac

if __name__ == "__main__":
    project = signac.get_project("../")
    for job in project:
        print(job.ws)
        job.document.pop("RDF_first_peak")
