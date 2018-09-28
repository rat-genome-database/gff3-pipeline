package edu.mcw.rgd.gff3;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.GeneQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Relation;
import edu.mcw.rgd.datamodel.ontologyx.TermWithStats;
import edu.mcw.rgd.process.Utils;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.object.MappingSqlQuery;

import javax.sql.DataSource;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * Date: 2/26/14
 * <p>
 * wrapper for DAO code
 */
public class RgdGff3Dao {

    private AliasDAO aliasDAO = new AliasDAO();
    private AnnotationDAO annotationDAO = new AnnotationDAO();
    private GeneDAO geneDAO = new GeneDAO();
    private QTLDAO qtlDAO = new QTLDAO();
    private MapDAO mapDAO = new MapDAO();
    private OntologyXDAO ontologyXDAO = new OntologyXDAO();
    private SequenceDAO sequenceDAO = new SequenceDAO();
    private SSLPDAO sslpDAO = new SSLPDAO();
    private StrainDAO strainDAO = new StrainDAO();
    private TranscriptDAO trDao = new TranscriptDAO();
    private VariantInfoDAO variantInfoDAO = new VariantInfoDAO();
    private XdbIdDAO xdbDao = new XdbIdDAO();


    /**
     * given rgd, id, return a gene object; for given gene, build description using the DescriptionGenerator
     * and set the gene.description to this value;
     * use cache to reduce number of database round trips
     *
     * @param geneRgdId gene rgd id
     * @return Gene object
     * @throws Exception
     */
    public Gene getGeneWithDescription(int geneRgdId) throws Exception {
        // first try to get object from cache
        Gene gene = _cacheGenes.get(geneRgdId);
        if (gene == null) {
            gene = geneDAO.getGene(geneRgdId);
            _cacheGenes.put(geneRgdId, gene);
            if (gene != null)
                gene.setDescription(Utils.getGeneDescription(gene));
        }
        return gene;
    }

    static private Map<Integer, Gene> _cacheGenes = new HashMap<>();

    public List<Gene> getActiveGenes(String chr, long startPos, long stopPos, int mapKey) throws Exception {
        return geneDAO.getActiveGenes(chr, startPos, stopPos, mapKey);
    }

    public List<Gene> getActiveGenes(int speciesTypeKey) throws Exception {
        return geneDAO.getActiveGenes(speciesTypeKey);
    }

    /**
     * get TermWithStats object given term acc; use term cache behind the scenes
     *
     * @param termAcc term acc id
     * @return TermWithStats object
     * @throws Exception
     */
    public TermWithStats getTerm(String termAcc) throws Exception {
        return ontologyXDAO.getTermWithStatsCached(termAcc);
    }


    List<Annotation> getAnnotations(int rgdId) throws Exception {
        return annotationDAO.getAnnotations(rgdId);
    }

    List<Annotation> getAnnotationsByAspect(int rgdId, String aspect) throws Exception {
        return annotationDAO.getAnnotationsByAspect(rgdId, aspect);
    }

    /**
     * get all direct descendant (child) terms of given term
     *
     * @param termAcc child term accession id
     * @return map of descendant terms: acc ids are the keys, term relations are the values
     * @throws Exception if something wrong happens in spring framework
     */
    public Map<String, Relation> getTermDescendants(String termAcc) throws Exception {
        return ontologyXDAO.getTermDescendants(termAcc);
    }

    List<SSLP> getActiveSslps(int speciesTypeKey) throws Exception {
        return sslpDAO.getActiveSSLPs(speciesTypeKey);
    }

    List<Sequence> getSslpSequences(int rgdId) throws Exception {
        return sequenceDAO.getObjectSequences(rgdId);
    }

    public Gene getGeneBySslpKey(int sslpKey) throws Exception {
        String sql = "select g.*,0 species_type_key from GENES g where GENE_KEY IN (SELECT GENE_KEY FROM RGD_GENE_SSLP WHERE SSLP_KEY=?)";
        GeneQuery q = new GeneQuery(sslpDAO.getDataSource(), sql);
        q.declareParameter(new SqlParameter(Types.INTEGER));
        q.compile();
        List<Gene> genes = q.execute(sslpKey);
        return genes == null || genes.isEmpty() ? null : genes.get(0);
    }

    public List<Alias> getAliases(int rgdId) throws Exception {
        return aliasDAO.getAliases(rgdId);
    }

    public List<MapData> getMapData(int rgdId, int mapKey) throws Exception {
        return mapDAO.getMapData(rgdId, mapKey);
    }

    /**
     * return all positions for given map, chromosome and object type
     *
     * @param chromosome chromosome
     * @param mapKey     map key
     * @param objectKey  object key
     * @return List of MapData objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<MapData> getMapDataByMapKeyChr(String chromosome, int mapKey, int objectKey) throws Exception {
        return mapDAO.getMapDataByMapKeyChr(chromosome, mapKey, objectKey);
    }

    public Map<String, Integer> getChromosomeSizes(int mapKey) throws Exception {
        return mapDAO.getChromosomeSizes(mapKey);
    }

    public List<Transcript> getTranscriptsForGene(int geneRgdId) throws Exception {
        return trDao.getTranscriptsForGene(geneRgdId);
    }

    public List<TranscriptFeature> getFeatures(int rgdId, int mapKey) throws Exception {
        return trDao.getFeatures(rgdId, mapKey);
    }

    public List<VariantInfo> getClinVarVariants() throws Exception {
        return variantInfoDAO.getVariantsBySource("CLINVAR");
    }

    public List<XdbId> getXdbIds(XdbId filter) throws Exception {
        return xdbDao.getXdbIds(filter);
    }

    public Strain getStrain(int rgdId) throws Exception {
        return strainDAO.getStrain(rgdId);
    }

    public List<Strain> getMappedStrains(int mapKey) throws Exception {
        return strainDAO.getMappedStrains(mapKey);
    }

    public List<QTL> getActiveQTLs(int speciesTypeKey) throws Exception {
        return qtlDAO.getActiveQTLs(speciesTypeKey);
    }

    public List<DbSnp> getDbSnps(int mapKey, String source) throws Exception {
        String sql = "SELECT * FROM db_snp WHERE map_key=? AND source=?";
        DbSnpQuery q = new DbSnpQuery(geneDAO.getDataSource(), sql);
        return geneDAO.execute(q, mapKey, source);
    }

    class DbSnp {
        public String snpName; // rsXXX
        public String source;  // dbSnpXXX
        public int mapKey;
        public String chr;
        public int pos;
        public String snpClass; // 'snp'
        public String allele;
        public String refAllele;
        public String functionClass;
        public String clinicalSignificance;
    }

    public class DbSnpQuery extends MappingSqlQuery {
        public DbSnpQuery(DataSource ds, String query) {
            super(ds, query);
        }

        protected Object mapRow(ResultSet rs, int rowNum) throws SQLException {
            DbSnp obj = new DbSnp();
            obj.snpName = rs.getString("snp_name");
            obj.source = rs.getString("source");
            obj.mapKey = rs.getInt("map_key");
            obj.chr = rs.getString("chromosome");
            obj.pos = rs.getInt("position");
            obj.snpClass = rs.getString("snp_class");
            obj.allele = rs.getString("allele");
            obj.refAllele = rs.getString("ref_allele");
            obj.functionClass = rs.getString("function_class");
            obj.clinicalSignificance = rs.getString("clinical_significance");

            return obj;
        }
    }
}