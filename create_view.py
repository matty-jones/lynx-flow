import signac
import os

if __name__ == "__main__":
    cwd = os.getcwd()
    project = signac.get_project("./")
    schema = project.detect_schema()
    for dimensions in schema["dimensions"].values()[0]:
        for z_size in schema["z_reactor_size"].values()[0]:
            job_list = project.find_job_ids(
                {
                    "dimensions": dimensions,
                    "z_reactor_size": z_size,
                    "job_type": "child",
                }
            )
            view_name = "".join(["layer", str(layer), "z_size", str(z_size)])
            project.create_linked_view(
                prefix=os.path.join(cwd, "views", view_name), job_ids=job_list
            )
