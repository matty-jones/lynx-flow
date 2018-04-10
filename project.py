import flow


class LynxProject(flow.FlowProject):
    def __init__(self, *args, **kwargs):
        super(LynxProject, self).__init__(*args, **kwargs)
        self.add_operation(
            name='generated',
            cmd=lambda job: "python operations.py generate {}".format(job),
            post=[LynxProject.generated]
        )

    @flow.staticlabel()
    def generated(job):
        return job.isfile('output.hoomdxml')

    @flow.staticlabel()
    def generated(job):
        return job.isfile('output.gsd')


if __name__ == '__main__':
    LynxProject().main()
