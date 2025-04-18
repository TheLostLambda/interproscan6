import java.io.File
import java.nio.file.*
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process DOWNLOAD {
    maxForks 4
    label    'local', 'ips6_container'

    input:
    tuple val(name), val(arcname), val(version)
    val iprscan_version
    val outdir

    output:
    tuple val(name), val(version), val("${outdir}/${arcname}/${version}")

    script:
    """
    cd ${outdir}
    curl -OJ ${InterProScan.FTP_URL}/${iprscan_version}/${arcname}/${arcname}-${version}.tar.gz
    curl -OJ ${InterProScan.FTP_URL}/${iprscan_version}/${arcname}/${arcname}-${version}.tar.gz.md5
    md5sum -c ${arcname}-${version}.tar.gz.md5 || { echo "Error: MD5 checksum failed" >&2; exit 1; }
    tar -zxf ${arcname}-${version}.tar.gz
    rm ${arcname}-${version}.tar.gz*
    chmod 777 -R ${arcname}
    """
}

process FIND_MISSING_APP_DATA {
    input:
    tuple val(n), val(v), val(p)  // state dependency
    val json_database
    val apps_to_run
    val app_dirs
    val datadir

    output:
    val with_data,          emit: with_data
    val without_data,       emit: without_data

    exec:
    File file = new File(json_database)
    def json = new JsonSlurper().parse(file)
    def normalised_json = [:]
    json.each { key, value ->
        normalised_json[key.replaceAll(/[\s\-]+/, '').toLowerCase()] = value
    }

    with_data = [] as Set
    without_data = [] as Set
    apps_to_run.each { db_name ->
        if (app_dirs.containsKey(db_name)) {
            def normalised_name = db_name.replaceAll(/[\s\-]+/, '').toLowerCase()
            assert normalised_json.containsKey(normalised_name)
            def db_version = normalised_json[normalised_name]
            def db_dir_parts = app_dirs[normalised_name]
            def db_dir = db_dir_parts[0]
            def db_subdir = db_dir_parts[1]
            Path path = Paths.get("${datadir}/${db_dir}/${db_version}/${db_subdir}")
            if (Files.exists(path)) {
                with_data.add( [ normalised_name, db_version, path.toString() ])
            } else {
                without_data.add( [ normalised_name, db_dir, db_version ] )
            }
        }
    }
}

process WAIT_FOR_DOWNLOADS {
    cache false  // Stops the esotericsoftware.kryo.serializers warning

    input:
    val list_databases

    output:
    val map_databases

    exec:
    map_databases = list_databases.collectEntries { name, version, dirpath ->
        [(name): [version: version, dirpath: dirpath]]
    } 
}
