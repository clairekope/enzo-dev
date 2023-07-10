import yt
es = yt.load_simulation('TestStarParticleSingle.enzo', 'Enzo')
es.get_time_series()
for ds in es:
    plt = yt.ProjectionPlot(ds, 'x', 'density')
    plt.annotate_particles(width=(400,"pc"))
    plt.save()
    yt.ProjectionPlot(ds, 'x', 'temperature', weight_field='density').save()
