import signac
import os

if __name__ == "__main__":
    cwd = os.getcwd()
    project = signac.get_project("./")
    for layer in [1, 2, 3]:
        for z_size in [10.0, 15.0, 20.0, 25.0]:
            job_list = project.find_job_ids(
                {
                    "dimensions": "".join(["10x10x", str(layer)]),
                    "z_reactor_size": z_size,
                }
            )
            view_name = "".join(["layer", str(layer), "z_size", str(z_size)])
            project.create_linked_view(
                prefix=os.path.join(cwd, "views", view_name), job_ids=job_list
            )
