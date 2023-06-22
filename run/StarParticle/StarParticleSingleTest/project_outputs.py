import yt
es = yt.load_simulation('TestStarParticleSingle.enzo', 'Enzo')
es.get_time_series()
for ds in es:
    yt.ProjectionPlot(ds, 'x', 'density').save()
    yt.ProjectionPlot(ds, 'x', 'temperature', weight_field='density').save()
